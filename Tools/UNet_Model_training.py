#!/usr/bin/env python3
"""
Fine-tune the NHP-BrainExtraction 2D U-Net (Wang et al., 2021, NeuroImage 235:118001)
on the latest EDNiX outputs (anat + func), starting from model-02-_cross_species-epoch.

Run:  python3 Tools/UNet_Model_training.py            (everything: discover -> split -> cache -> train -> test)
      python3 Tools/UNet_Model_training.py --dry-run  (discovery + split only, prints counts)

Training pairs (current EDNiX naming, native space, same voxel grid):
  anat : <ses>/anat/native/volumes/{ID}_space-acpc_desc-template_{T1w|T2w}.nii.gz      (whole head, N4)
         <ses>/anat/native/volumes/masks/{ID}_backtonative_mask.nii.gz                (manual, priority)
      or <ses>/anat/native/volumes/masks/{ID}_space-acpc_mask.nii.gz                  (final automatic)
  func : <ses>/func/raw/preprocessing/all_runs_space-func_desc-fMRI_Mean_Image.nii.gz
         <ses>/func/raw/masks/{ID}_fMRI_mask.nii.gz  (= {ID}_final_mask resampled when a manual mask exists)

The checkpoint written (checkpoints/model-cross_species_ft-best-epoch) keeps the original
architecture and checkpoint layout: it is a drop-in replacement for muSkullStrip.py -model.
"""
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import os
import random
import sys
import time
import types
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from pathlib import Path
import multiprocessing as mp

import numpy as np
import nibabel as nib
import scipy.ndimage as snd

# ----------------------------------------------------------------------------------------
# Default locations (first existing candidate wins; every one can be overridden on the CLI)
# ----------------------------------------------------------------------------------------
EDNIX_ROOT = Path(__file__).resolve().parents[1]
UNET_DIR_CANDIDATES = [
    EDNIX_ROOT / 'Singularity_library' / 'Singularity' / 'NHP-BrainExtraction' / 'UNet_Model',
    EDNIX_ROOT / 'Tool_library' / 'Singularity' / 'NHP-BrainExtraction' / 'UNet_Model',
    Path('/home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/NHP-BrainExtraction/UNet_Model'),
    Path('/srv/projects/easymribrain/code/EDNiX/Singularity_library/Singularity/NHP-BrainExtraction/UNet_Model'),
    Path('/srv/projects/easymribrain/code/EDNiX/Tool_library/Singularity/NHP-BrainExtraction/UNet_Model'),
]
DATA_ROOT_CANDIDATES = [Path('/scratch2/EDNiX'), Path('/srv/projects/easymribrain/data/MRI')]
OUT_DIR_CANDIDATES = [
    (Path('/srv/projects/easymribrain/scratch'), 'Training_UNet_cross_species'),
    (Path('/scratch2'), 'Training_UNet_cross_species'),
]
INIT_MODEL_NAME = 'model-02-_cross_species-epoch'

RESCALE_DIM = 256   # must stay = muSkullStrip default (-rescale 256)
NUM_SLICE = 3       # must stay = muSkullStrip default (-slice 3)
CONV_BLOCK = 5
KERNEL_ROOT = 16
CACHE_VERSION = 'v1'


def first_existing(paths):
    for p in paths:
        if Path(p).exists():
            return Path(p)
    return None


def default_out_dir():
    for parent, name in OUT_DIR_CANDIDATES:
        if parent.is_dir() and os.access(parent, os.W_OK):
            return parent / name
    return EDNIX_ROOT.parent / 'Training_UNet_cross_species'


# ----------------------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------------------
_LOG_FH = None


def log(msg):
    line = time.strftime('%H:%M:%S ') + str(msg)
    print(line, flush=True)
    if _LOG_FH is not None:
        _LOG_FH.write(line + '\n')
        _LOG_FH.flush()


# ----------------------------------------------------------------------------------------
# 1. Discovery of image/mask pairs
# ----------------------------------------------------------------------------------------
def find_subject_dirs(roots, max_depth=4):
    seen = set()
    for root in roots:
        for depth in range(max_depth + 1):
            pattern = os.path.join(str(root), *(['*'] * depth), 'sub-*')
            for d in glob(pattern):
                if not os.path.isdir(d):
                    continue
                parts = Path(d).parts
                if any(p.lower() in ('archives', 'archive', 'old', 'derivatives_old') for p in parts):
                    continue
                real = os.path.realpath(d)
                if real not in seen:
                    seen.add(real)
                    yield Path(d)


def species_dataset(sub_dir):
    parts = sub_dir.parts
    for i in range(len(parts) - 1, 0, -1):
        if parts[i].upper().startswith('BIDS'):
            return parts[i - 1], parts[i]
    return sub_dir.parent.parent.name, sub_dir.parent.name


def same_grid(img_path, msk_path):
    try:
        a, b = nib.load(str(img_path)), nib.load(str(msk_path))
    except Exception as e:  # corrupted file
        return False, f'unreadable ({e})'
    sa, sb = a.shape[:3], b.shape[:3]
    if sa != sb:
        return False, f'shape {sa} vs {sb}'
    if len(a.shape) > 3 and int(np.prod(a.shape[3:])) > 1:
        return False, f'image is 4D {a.shape}'
    if not np.allclose(a.affine, b.affine, atol=1e-2):
        return False, 'affine mismatch'
    return True, ''


def discover_pairs(roots, func_masks='all', anat_contrasts=('T1w', 'T2w')):
    pairs, skipped = [], []
    for sub_dir in find_subject_dirs(roots):
        species, dataset = species_dataset(sub_dir)
        sub = sub_dir.name
        ses_dirs = sorted(Path(p) for p in glob(str(sub_dir / 'ses-*')) if os.path.isdir(p)) or [sub_dir]
        for ses_dir in ses_dirs:
            ses = ses_dir.name if ses_dir != sub_dir else 'ses-NA'

            # ---- anat -------------------------------------------------------------
            vol_dir = ses_dir / 'anat' / 'native' / 'volumes'
            msk_dir = vol_dir / 'masks'
            ids = set()
            for f in glob(str(msk_dir / '*_space-acpc_mask.nii.gz')):
                ids.add(os.path.basename(f)[:-len('_space-acpc_mask.nii.gz')])
            for f in glob(str(msk_dir / '*_backtonative_mask.nii.gz')):
                ids.add(os.path.basename(f)[:-len('_backtonative_mask.nii.gz')])
            for ID in sorted(ids):
                manual = msk_dir / f'{ID}_backtonative_mask.nii.gz'
                auto = msk_dir / f'{ID}_space-acpc_mask.nii.gz'
                mask = manual if manual.is_file() else auto
                for T in anat_contrasts:
                    img = vol_dir / f'{ID}_space-acpc_desc-template_{T}.nii.gz'
                    if not img.is_file():
                        continue
                    ok, why = same_grid(img, mask)
                    rec = dict(species=species, dataset=dataset, subject=sub, session=ses, kind='anat',
                               contrast=T, manual_mask=int(mask == manual), img=str(img), mask=str(mask),
                               mtime=max(img.stat().st_mtime, mask.stat().st_mtime))
                    (pairs if ok else skipped).append(rec if ok else {**rec, 'reason': why})

            # ---- func -------------------------------------------------------------
            raw = ses_dir / 'func' / 'raw'
            mean_img = raw / 'preprocessing' / 'all_runs_space-func_desc-fMRI_Mean_Image.nii.gz'
            if mean_img.is_file():
                for f in sorted(glob(str(raw / 'masks' / '*_fMRI_mask.nii.gz'))):
                    ID = os.path.basename(f)[:-len('_fMRI_mask.nii.gz')]
                    manual = int((raw / 'masks' / f'{ID}_final_mask.nii.gz').is_file())
                    if func_masks == 'manual' and not manual:
                        continue
                    ok, why = same_grid(mean_img, f)
                    rec = dict(species=species, dataset=dataset, subject=sub, session=ses, kind='func',
                               contrast='EPI', manual_mask=manual, img=str(mean_img), mask=f,
                               mtime=max(mean_img.stat().st_mtime, os.path.getmtime(f)))
                    (pairs if ok else skipped).append(rec if ok else {**rec, 'reason': why})

    # Same data mirrored under several roots: keep the most recent outputs only.
    best = {}
    for r in pairs:
        key = (r['species'].lower(), r['dataset'], r['subject'], r['session'], r['kind'], r['contrast'])
        if key not in best or r['mtime'] > best[key]['mtime']:
            best[key] = r
    pairs = sorted(best.values(), key=lambda r: (r['species'], r['dataset'], r['subject'], r['session'],
                                                 r['kind'], r['contrast']))
    for r in pairs:
        r['subject_key'] = f"{r['species']}/{r['dataset']}/{r['subject']}"
        r['stratum'] = f"{r['dataset']}|{r['kind']}"
        r['uid'] = hashlib.blake2b(f"{r['img']}|{r['mask']}".encode(), digest_size=8).hexdigest()
    return pairs, skipped


# ----------------------------------------------------------------------------------------
# 2. Subject-level split, stratified by dataset
# ----------------------------------------------------------------------------------------
def make_split(pairs, val_ratio, test_ratio, seed):
    rng = np.random.default_rng(seed)
    by_ds = defaultdict(set)
    for r in pairs:
        by_ds[(r['species'], r['dataset'])].add(r['subject_key'])
    split = {}
    for ds in sorted(by_ds):
        subs = sorted(by_ds[ds])
        subs = [subs[i] for i in rng.permutation(len(subs))]
        n = len(subs)
        n_val = int(round(n * val_ratio))
        n_test = int(round(n * test_ratio))
        if n >= 3 and val_ratio > 0 and n_val == 0:
            n_val = 1
        if n >= 4 and test_ratio > 0 and n_test == 0:
            n_test = 1
        while n - n_val - n_test < 1:
            if n_test > 0:
                n_test -= 1
            elif n_val > 0:
                n_val -= 1
        for i, s in enumerate(subs):
            split[s] = 'test' if i < n_test else ('val' if i < n_test + n_val else 'train')
    return split


def load_or_make_split(pairs, path, val_ratio, test_ratio, seed, resplit):
    if path.is_file() and not resplit:
        split = json.loads(path.read_text())
        new = sorted({r['subject_key'] for r in pairs} - set(split))
        if new:
            log(f'WARNING: {len(new)} subjects absent from {path.name} -> assigned to train '
                f'(use --resplit to redo a stratified split).')
            for s in new:
                split[s] = 'train'
    else:
        split = make_split(pairs, val_ratio, test_ratio, seed)
    path.write_text(json.dumps(split, indent=1, sort_keys=True))
    for r in pairs:
        r['split'] = split[r['subject_key']]
    return split


# ----------------------------------------------------------------------------------------
# 3. Cache: volumes normalised/rescaled exactly like dataset.VolumeDataset + BlockDataset
# ----------------------------------------------------------------------------------------
def _cache_one(args):
    import torch
    import torch.nn.functional as F
    torch.set_num_threads(1)
    rec, cache_root, min_frac, max_frac = args
    d = Path(cache_root) / rec['uid']
    meta_f = d / 'meta.json'
    key = f"{CACHE_VERSION}|{rec['img']}|{rec['mask']}|{rec['mtime']}"
    if meta_f.is_file():
        meta = json.loads(meta_f.read_text())
        if meta.get('key') == key:
            return rec['uid'], meta, ''
    try:
        img = np.asarray(nib.load(rec['img']).get_fdata(dtype=np.float32)).squeeze()
        msk = np.asarray(nib.load(rec['mask']).dataobj).squeeze() > 0
    except Exception as e:
        return rec['uid'], None, f'read error: {e}'
    if img.ndim != 3:
        return rec['uid'], None, f'ndim={img.ndim}'
    img[~np.isfinite(img)] = 0
    lo, hi = float(img.min()), float(img.max())
    if hi <= lo:
        return rec['uid'], None, 'constant image'
    frac = float(msk.mean())
    if not (min_frac <= frac <= max_frac):
        return rec['uid'], None, f'mask fraction {frac:.4f} out of [{min_frac}, {max_frac}]'
    img = (img - lo) / (hi - lo)                                   # same 0-1 min-max as VolumeDataset
    raw_shape = list(img.shape)
    factor = float(RESCALE_DIM) / float(max(raw_shape))            # same factor as BlockDataset
    with torch.no_grad():
        ti = F.interpolate(torch.from_numpy(img)[None, None], scale_factor=factor,
                           mode='trilinear', align_corners=False)[0, 0]
        tm = F.interpolate(torch.from_numpy(msk.astype(np.float32))[None, None], scale_factor=factor,
                           mode='nearest')[0, 0]
    ri, rm = ti.numpy().astype(np.float16), (tm.numpy() > 0.5).astype(np.uint8)
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / 'img.npy', ri)
    np.save(d / 'msk.npy', rm)
    # starts s of 3-slice blocks whose centre slice (s+1) contains brain, per axis
    fg = []
    for ax in range(3):
        other = tuple(a for a in range(3) if a != ax)
        has = rm.any(axis=other)
        fg.append([int(s) for s in range(rm.shape[ax] - NUM_SLICE + 1) if has[s + 1]])
    meta = dict(key=key, raw_shape=raw_shape, rescale_shape=list(ri.shape), factor=factor,
                mask_fraction=frac, fg_starts=fg)
    meta_f.write_text(json.dumps(meta))
    return rec['uid'], meta, ''


def build_cache(pairs, cache_root, n_workers, min_frac, max_frac):
    cache_root.mkdir(parents=True, exist_ok=True)
    jobs = [(r, str(cache_root), min_frac, max_frac) for r in pairs]
    metas, failed = {}, []
    t0 = time.time()
    ctx = mp.get_context('spawn')
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        for i, (uid, meta, err) in enumerate(ex.map(_cache_one, jobs, chunksize=1)):
            if meta is None:
                failed.append((uid, err))
            else:
                metas[uid] = meta
            if (i + 1) % 50 == 0 or i + 1 == len(jobs):
                log(f'  cache {i + 1}/{len(jobs)}  ({time.time() - t0:.0f}s)')
    return metas, failed


# ----------------------------------------------------------------------------------------
# 4. Model helpers
# ----------------------------------------------------------------------------------------
def load_unet_module(unet_dir):
    spec = importlib.util.spec_from_file_location('nhp_unet_model', str(unet_dir / 'model.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _per_sample_bn_forward(self, x):
    import torch.nn.functional as F
    # = BatchNorm2d in training mode with batch size 1, which is what muSkullStrip.py actually
    # runs at inference (it never calls model.eval() and feeds one block at a time).
    return F.instance_norm(x, None, None, self.weight, self.bias, True, 0.0, self.eps)


def patch_bn_per_sample(model):
    import torch.nn as nn
    n = 0
    for m in model.modules():
        if isinstance(m, nn.BatchNorm2d):
            m.forward = types.MethodType(_per_sample_bn_forward, m)
            n += 1
    return n


def torch_load(path):
    import torch
    try:
        return torch.load(path, map_location='cpu', weights_only=False)
    except TypeError:  # torch < 1.13
        return torch.load(path, map_location='cpu')


def build_model(unet_mod, init_path, device):
    model = unet_mod.UNet2d(dim_in=NUM_SLICE, num_conv_block=CONV_BLOCK, kernel_root=KERNEL_ROOT)
    ckpt = torch_load(str(init_path))
    state = ckpt['state_dict'] if isinstance(ckpt, dict) and 'state_dict' in ckpt else ckpt
    state = OrderedDict((k.replace('module.', '', 1), v) for k, v in state.items())
    model.load_state_dict(state)
    patch_bn_per_sample(model)
    return model.to(device)


# ----------------------------------------------------------------------------------------
# 5. Training data: volume-grouped random 3-slice blocks
# ----------------------------------------------------------------------------------------
def make_block(vol, ax, s):
    """Same slicing/permutation/padding as dataset.BlockDataset.__getitem__ (numpy)."""
    if ax == 0:
        b = vol[s:s + NUM_SLICE, :, :]
    elif ax == 1:
        b = np.transpose(vol[:, s:s + NUM_SLICE, :], (1, 0, 2))
    else:
        b = np.transpose(vol[:, :, s:s + NUM_SLICE], (2, 0, 1))
    return b


class BlockVisitDataset:
    """Item = (volume index, [(axis, start), ...]) -> stacked blocks from one volume."""

    def __init__(self, cache_root, uids, gamma_p, gamma_range, lru=48):
        self.cache_root, self.uids = str(cache_root), uids
        self.gamma_p, self.gamma_range = gamma_p, gamma_range
        self.lru, self._open = lru, OrderedDict()

    def _arrays(self, uid):
        if uid in self._open:
            self._open.move_to_end(uid)
            return self._open[uid]
        d = os.path.join(self.cache_root, uid)
        arr = (np.load(os.path.join(d, 'img.npy'), mmap_mode='r'),
               np.load(os.path.join(d, 'msk.npy'), mmap_mode='r'))
        self._open[uid] = arr
        if len(self._open) > self.lru:
            self._open.popitem(last=False)
        return arr

    def __len__(self):
        return len(self.uids)

    def __getitem__(self, item):
        import torch
        v, draws, seed = item
        rng = np.random.default_rng(seed)
        img, msk = self._arrays(self.uids[v])
        X = np.zeros((len(draws), NUM_SLICE, RESCALE_DIM, RESCALE_DIM), np.float32)
        Y = np.zeros((len(draws), RESCALE_DIM, RESCALE_DIM), np.int64)
        for i, (ax, s) in enumerate(draws):
            bi = make_block(img, ax, s).astype(np.float32)
            bm = make_block(msk, ax, s)[1]
            if rng.random() < self.gamma_p:
                bi = np.power(np.clip(bi, 0, 1), rng.uniform(*self.gamma_range))
            h, w = bi.shape[1:]
            X[i, :, :h, :w] = bi
            Y[i, :h, :w] = bm
        return torch.from_numpy(X), torch.from_numpy(Y)


def collate_visits(batch):
    import torch
    return torch.cat([b[0] for b in batch]), torch.cat([b[1] for b in batch])


class EpochVisitSampler:
    def __init__(self, train_idx, weights, metas_list, blocks_per_volume, blocks_per_visit, fg_p, seed):
        self.train_idx = np.asarray(train_idx)
        w = np.asarray(weights, float)
        self.p = w / w.sum()
        self.metas = metas_list
        self.n_visits = max(1, int(round(len(train_idx) * blocks_per_volume / blocks_per_visit)))
        self.k, self.fg_p = blocks_per_visit, fg_p
        self.rng = np.random.default_rng(seed)

    def __len__(self):
        return self.n_visits

    def __iter__(self):
        vols = self.rng.choice(self.train_idx, size=self.n_visits, replace=True, p=self.p)
        for v in vols:
            m = self.metas[v]
            draws = []
            for _ in range(self.k):
                ax = int(self.rng.integers(3))
                n_starts = m['rescale_shape'][ax] - NUM_SLICE + 1
                fg = m['fg_starts'][ax]
                if fg and self.rng.random() < self.fg_p:
                    s = int(fg[self.rng.integers(len(fg))])
                else:
                    s = int(self.rng.integers(n_starts))
                draws.append((ax, s))
            yield int(v), draws, int(self.rng.integers(2 ** 31))


# ----------------------------------------------------------------------------------------
# 6. Evaluation = muSkullStrip/function.predict_volumes, batched on GPU
# ----------------------------------------------------------------------------------------
def extract_large_comp(m):
    labs, n = snd.label(m)
    if n == 0:
        return m.astype(bool)
    c = np.bincount(labs.ravel())
    c[0] = 0
    return labs == c.argmax()


def fill_holes(m):
    return ~extract_large_comp(~m)


def predict_volume(model, img_resc, raw_shape, device, batch, amp_dtype):
    import torch
    import torch.nn.functional as F
    vol = torch.from_numpy(np.asarray(img_resc, dtype=np.float32)).to(device)
    R = vol.shape
    acc = torch.zeros(R, device=device)
    perms = {0: ((0, 1, 2), (0, 1, 2)), 1: ((1, 0, 2), (1, 0, 2)), 2: ((2, 0, 1), (1, 2, 0))}
    with torch.no_grad():
        for ax in range(3):
            fwd, back = perms[ax]
            V = vol.permute(*fwd)
            D, H, W = V.shape
            blocks = V.unfold(0, NUM_SLICE, 1).permute(0, 3, 1, 2)  # (D-2, 3, H, W)
            blocks = F.pad(blocks, (0, RESCALE_DIM - W, 0, RESCALE_DIM - H))
            pr = torch.zeros((D, H, W), device=device)  # first/last slice stay 0, as in predict_volumes
            for i in range(0, blocks.shape[0], batch):
                xb = blocks[i:i + batch].contiguous()
                with torch.autocast(device_type=device.type, dtype=amp_dtype, enabled=amp_dtype is not None):
                    out = model(xb)
                p = torch.softmax(out.float(), dim=1)[:, 1, :H, :W]
                pr[i + 1:i + 1 + p.shape[0]] = p
            acc += pr.permute(*back)
        prob = acc / 3.0
        prob = F.interpolate(prob[None, None], size=tuple(raw_shape), mode='trilinear',
                             align_corners=False)[0, 0]
    m = (prob > 0.5).cpu().numpy()
    return fill_holes(extract_large_comp(m))


def dice(a, b):
    s = a.sum() + b.sum()
    return float(2.0 * np.logical_and(a, b).sum() / s) if s > 0 else 1.0


class Evaluator:
    def __init__(self, recs, cache_root, device, batch, amp_dtype):
        self.recs, self.cache_root = recs, Path(cache_root)
        self.device, self.batch, self.amp_dtype = device, batch, amp_dtype
        self._gt = {}

    def _gt_mask(self, r):
        if r['uid'] not in self._gt:
            m = np.asarray(nib.load(r['mask']).dataobj).squeeze() > 0
            self._gt[r['uid']] = (np.packbits(m.ravel()), m.shape)
        bits, shape = self._gt[r['uid']]
        return np.unpackbits(bits, count=int(np.prod(shape))).reshape(shape).astype(bool)

    def run(self, model):
        model.eval()
        out = {}
        for r in self.recs:
            img = np.load(self.cache_root / r['uid'] / 'img.npy', mmap_mode='r')
            gt = self._gt_mask(r)
            pr = predict_volume(model, img, gt.shape, self.device, self.batch, self.amp_dtype)
            out[r['uid']] = dice(gt, pr)
        model.train()
        return out


def macro_dice(recs, d):
    by = defaultdict(list)
    for r in recs:
        if r['uid'] in d:
            by[r['stratum']].append(d[r['uid']])
    per = {k: float(np.mean(v)) for k, v in sorted(by.items())}
    return (float(np.mean(list(per.values()))) if per else float('nan')), per


# ----------------------------------------------------------------------------------------
# 7. Main
# ----------------------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-roots', nargs='+', default=None,
                   help='EDNiX data roots (default: existing among %s)' % [str(x) for x in DATA_ROOT_CANDIDATES])
    p.add_argument('--unet-dir', default=None, help='NHP-BrainExtraction/UNet_Model directory')
    p.add_argument('--init-model', default=None, help=f'init checkpoint (default: {INIT_MODEL_NAME} in UNet_Model/models)')
    p.add_argument('--out-dir', default=None, help='output dir (default: %s)' % default_out_dir())
    p.add_argument('--func-masks', choices=['all', 'manual'], default='all',
                   help="'manual': only func sessions with a hand-corrected {ID}_final_mask")
    p.add_argument('--anat-contrasts', nargs='+', default=['T1w', 'T2w'])
    p.add_argument('--val-ratio', type=float, default=0.15)
    p.add_argument('--test-ratio', type=float, default=0.15)
    p.add_argument('--resplit', action='store_true', help='ignore an existing split.json')
    p.add_argument('--seed', type=int, default=1234)
    p.add_argument('--epochs', type=int, default=30)
    p.add_argument('--patience', type=int, default=6, help='early stopping on validation macro-Dice')
    p.add_argument('--lr', type=float, default=1e-4)
    p.add_argument('--min-lr', type=float, default=1e-6)
    p.add_argument('--blocks-per-volume', type=int, default=48,
                   help='3-slice blocks drawn per training volume and epoch (original script: all, ~750)')
    p.add_argument('--blocks-per-visit', type=int, default=8)
    p.add_argument('--batch-size', type=int, default=32)
    p.add_argument('--fg-p', type=float, default=0.8, help='prob. of drawing a block centred on brain')
    p.add_argument('--balance-power', type=float, default=0.5,
                   help='stratum (dataset x anat/func) sampling mass ∝ n^power (1=natural, 0=fully balanced)')
    p.add_argument('--dice-weight', type=float, default=1.0, help='loss = CE + w * soft Dice')
    p.add_argument('--gamma-p', type=float, default=0.3)
    p.add_argument('--val-every', type=int, default=1)
    p.add_argument('--eval-batch', type=int, default=64)
    p.add_argument('--min-mask-frac', type=float, default=0.001)
    p.add_argument('--max-mask-frac', type=float, default=0.8)
    p.add_argument('--n-cpus', type=int, default=min(16, os.cpu_count() or 1))
    p.add_argument('--cpu', action='store_true', help='force CPU even if CUDA is available')
    p.add_argument('--dry-run', action='store_true', help='discovery + split only')
    return p.parse_args()


def main():
    global _LOG_FH
    args = parse_args()
    for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ.setdefault(v, str(args.n_cpus))

    roots = [Path(r) for r in args.data_roots] if args.data_roots else [p for p in DATA_ROOT_CANDIDATES if p.is_dir()]
    if not roots:
        sys.exit(f'No data root found among {DATA_ROOT_CANDIDATES}; pass --data-roots.')
    unet_dir = Path(args.unet_dir) if args.unet_dir else first_existing(UNET_DIR_CANDIDATES)
    if unet_dir is None or not (unet_dir / 'model.py').is_file():
        sys.exit(f'NHP-BrainExtraction/UNet_Model not found among {UNET_DIR_CANDIDATES}; pass --unet-dir.')
    if args.init_model:
        init_path = Path(args.init_model)
    else:
        cands = sorted(glob(str(unet_dir / 'models' / (INIT_MODEL_NAME + '*'))))
        init_path = Path(cands[0]) if cands else unet_dir / 'models' / INIT_MODEL_NAME
    if not init_path.is_file():
        sys.exit(f'Init model not found: {init_path}')

    out = Path(args.out_dir) if args.out_dir else default_out_dir()
    (out / 'checkpoints').mkdir(parents=True, exist_ok=True)
    _LOG_FH = open(out / 'train.log', 'a')
    log('=' * 90)
    log(f'data roots : {[str(r) for r in roots]}')
    log(f'UNet dir   : {unet_dir}')
    log(f'init model : {init_path}')
    log(f'output dir : {out}')

    # ---- discovery -----------------------------------------------------------------------
    t0 = time.time()
    pairs, skipped = discover_pairs(roots, args.func_masks, tuple(args.anat_contrasts))
    log(f'discovery: {len(pairs)} pairs, {len(skipped)} skipped ({time.time() - t0:.0f}s)')
    for s in skipped[:20]:
        log(f"  skip {s['kind']} {s['img']} : {s['reason']}")
    if len(skipped) > 20:
        log(f'  ... {len(skipped) - 20} more skipped (see skipped.csv)')
    if not pairs:
        sys.exit('No training pair found.')
    split = load_or_make_split(pairs, out / 'split.json', args.val_ratio, args.test_ratio, args.seed, args.resplit)

    fields = ['uid', 'split', 'species', 'dataset', 'subject', 'session', 'kind', 'contrast', 'manual_mask',
              'img', 'mask']
    with open(out / 'manifest.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(pairs)
    if skipped:
        with open(out / 'skipped.csv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['kind', 'img', 'mask', 'reason'], extrasaction='ignore')
            w.writeheader()
            w.writerows(skipped)

    tab = defaultdict(lambda: defaultdict(int))
    for r in pairs:
        tab[(r['species'], r['dataset'], r['kind'])][r['split']] += 1
    log(f"{'species':<12}{'dataset':<22}{'kind':<6}{'train':>7}{'val':>6}{'test':>6}{'manual':>8}")
    for k in sorted(tab):
        man = sum(r['manual_mask'] for r in pairs if (r['species'], r['dataset'], r['kind']) == k)
        log(f"{k[0]:<12}{k[1]:<22}{k[2]:<6}{tab[k]['train']:>7}{tab[k]['val']:>6}{tab[k]['test']:>6}{man:>8}")
    n_subj = defaultdict(int)
    for s, sp in split.items():
        n_subj[sp] += 1
    log(f"subjects: train {n_subj['train']}, val {n_subj['val']}, test {n_subj['test']}")
    if args.dry_run:
        log('dry run: stop.')
        return

    # ---- cache ---------------------------------------------------------------------------
    log('building cache (normalise + rescale once, float16) ...')
    metas, failed = build_cache(pairs, out / 'cache', args.n_cpus, args.min_mask_frac, args.max_mask_frac)
    for uid, err in failed:
        r = next(x for x in pairs if x['uid'] == uid)
        log(f"  drop {r['kind']} {r['img']} : {err}")
    pairs = [r for r in pairs if r['uid'] in metas]
    size_gb = sum(np.prod(metas[r['uid']]['rescale_shape']) * 3 for r in pairs) / 1e9
    log(f'cache ready: {len(pairs)} volumes, ~{size_gb:.1f} GB')

    train = [r for r in pairs if r['split'] == 'train']
    val = [r for r in pairs if r['split'] == 'val']
    test = [r for r in pairs if r['split'] == 'test']
    if not train:
        sys.exit('Empty training set.')

    # ---- torch setup ---------------------------------------------------------------------
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)
    torch.set_num_threads(args.n_cpus)
    device = torch.device('cuda' if torch.cuda.is_available() and not args.cpu else 'cpu')
    amp_dtype = None
    if device.type == 'cuda':
        torch.backends.cudnn.benchmark = True
        amp_dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16
    log(f'device: {device}' + (f' ({torch.cuda.get_device_name(0)}, AMP {amp_dtype})' if device.type == 'cuda' else
                                f' ({args.n_cpus} threads) — expect slow training, a GPU is strongly advised'))

    unet_mod = load_unet_module(unet_dir)
    init_model = build_model(unet_mod, init_path, device)
    model = copy.deepcopy(init_model)
    patch_bn_per_sample(model)  # deepcopy does not carry the bound-method patch reliably
    n_par = sum(p.numel() for p in model.parameters())
    log(f'UNet2d: {n_par / 1e6:.2f} M parameters')

    uids = [r['uid'] for r in pairs]
    idx_of = {u: i for i, u in enumerate(uids)}
    metas_list = [metas[u] for u in uids]
    train_idx = [idx_of[r['uid']] for r in train]
    n_str = defaultdict(int)
    for r in train:
        n_str[r['stratum']] += 1
    weights = [n_str[r['stratum']] ** (args.balance_power - 1.0) for r in train]
    log('sampling mass per stratum: ' + ', '.join(
        f"{k}={n * n ** (args.balance_power - 1.0) / sum(m * m ** (args.balance_power - 1.0) for m in n_str.values()):.3f}"
        for k, n in sorted(n_str.items())))

    ds = BlockVisitDataset(out / 'cache', uids, args.gamma_p, (0.7, 1.5))
    sampler = EpochVisitSampler(train_idx, weights, metas_list, args.blocks_per_volume, args.blocks_per_visit,
                                args.fg_p, args.seed)
    n_workers = max(0, min(8, args.n_cpus - 1))
    loader = DataLoader(ds, batch_size=max(1, args.batch_size // args.blocks_per_visit), sampler=sampler,
                        collate_fn=collate_visits, num_workers=n_workers, pin_memory=device.type == 'cuda',
                        persistent_workers=n_workers > 0, prefetch_factor=4 if n_workers > 0 else None)

    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs, eta_min=args.min_lr)
    use_scaler = amp_dtype == torch.float16
    scaler = torch.cuda.amp.GradScaler(enabled=use_scaler) if device.type == 'cuda' else None
    ce = nn.CrossEntropyLoss()

    val_eval = Evaluator(val, out / 'cache', device, args.eval_batch, amp_dtype)
    test_eval = Evaluator(test, out / 'cache', device, args.eval_batch, amp_dtype)

    def save_ckpt(path, epoch, vdice):
        torch.save({'epoch': int(epoch),
                    'state_dict': OrderedDict((k, v.detach().cpu()) for k, v in model.state_dict().items()),
                    'validate_dice': [float(x) for x in vdice],
                    'init_model': str(init_path),
                    'config': {k: (v if isinstance(v, (int, float, str, bool, type(None))) else str(v))
                               for k, v in vars(args).items()}}, str(path))

    best, best_epoch, bad = -1.0, -1, 0
    if val:
        t = time.time()
        d0 = val_eval.run(init_model)
        best, per = macro_dice(val, d0)
        log(f'epoch -1 (init model) val macro-Dice {best:.4f}  ({time.time() - t:.0f}s)')
        for k, v in per.items():
            log(f'    {k:<34}{v:.4f}')
        save_ckpt(out / 'checkpoints' / 'model-cross_species_ft-best-epoch', -1, list(d0.values()))
    else:
        log('WARNING: no validation subject -> no model selection/early stopping, last epoch is kept.')

    log_f = out / 'train_log.csv'
    new_log = not log_f.is_file()
    with open(log_f, 'a', newline='') as fh:
        lw = csv.writer(fh)
        if new_log:
            lw.writerow(['epoch', 'lr', 'loss', 'blocks', 'seconds', 'val_macro_dice'])
        for epoch in range(args.epochs):
            model.train()
            t, losses, nb = time.time(), [], 0
            for it, (x, y) in enumerate(loader):
                x = x.to(device, non_blocking=True)
                y = y.to(device, non_blocking=True)
                with torch.autocast(device_type=device.type, dtype=amp_dtype or torch.float32,
                                    enabled=amp_dtype is not None):
                    logits = model(x)
                logits = logits.float()
                loss = ce(logits, y)
                if args.dice_weight > 0:
                    p1 = torch.softmax(logits, 1)[:, 1]
                    yf = y.float()
                    sd = (2 * (p1 * yf).sum() + 1.0) / (p1.sum() + yf.sum() + 1.0)
                    loss = loss + args.dice_weight * (1.0 - sd)
                opt.zero_grad(set_to_none=True)
                if use_scaler:
                    scaler.scale(loss).backward()
                    scaler.step(opt)
                    scaler.update()
                else:
                    loss.backward()
                    opt.step()
                losses.append(float(loss.detach()))
                nb += x.shape[0]
                if (it + 1) % 100 == 0:
                    log(f'  epoch {epoch} it {it + 1}/{len(loader)} loss {np.mean(losses[-100:]):.4f} '
                        f'({nb / (time.time() - t):.0f} blocks/s)')
            lr = opt.param_groups[0]['lr']
            sched.step()
            dt = time.time() - t
            vm = float('nan')
            if val and ((epoch + 1) % args.val_every == 0 or epoch + 1 == args.epochs):
                tv = time.time()
                dv = val_eval.run(model)
                vm, per = macro_dice(val, dv)
                log(f'epoch {epoch} loss {np.mean(losses):.4f} lr {lr:.2e} train {dt:.0f}s | '
                    f'val macro-Dice {vm:.4f} ({time.time() - tv:.0f}s)')
                for k, v in per.items():
                    log(f'    {k:<34}{v:.4f}')
                if vm > best:
                    best, best_epoch, bad = vm, epoch, 0
                    save_ckpt(out / 'checkpoints' / 'model-cross_species_ft-best-epoch', epoch, list(dv.values()))
                    log(f'    -> new best model (epoch {epoch})')
                else:
                    bad += 1
            else:
                log(f'epoch {epoch} loss {np.mean(losses):.4f} lr {lr:.2e} train {dt:.0f}s')
            save_ckpt(out / 'checkpoints' / 'model-last-epoch', epoch, [])
            lw.writerow([epoch, lr, float(np.mean(losses)), nb, round(dt, 1), vm])
            fh.flush()
            if val and bad >= args.patience:
                log(f'early stopping: no improvement for {bad} validations.')
                break

    if not val:
        save_ckpt(out / 'checkpoints' / 'model-cross_species_ft-best-epoch', epoch, [])
    best_path = out / 'checkpoints' / 'model-cross_species_ft-best-epoch'
    if val and best_epoch < 0:
        log('WARNING: fine-tuning never beat the init model on validation; the "best" checkpoint is the init model.')

    # ---- final evaluation: init vs fine-tuned, on val and held-out test ---------------------
    final = build_model(unet_mod, best_path, device)
    rows = []
    for name, ev, recs in (('val', val_eval, val), ('test', test_eval, test)):
        if not recs:
            continue
        di, df = ev.run(init_model), ev.run(final)
        for r in recs:
            rows.append({**{k: r[k] for k in ('uid', 'species', 'dataset', 'subject', 'session', 'kind',
                                              'contrast', 'manual_mask')},
                         'split': name, 'dice_init': di[r['uid']], 'dice_finetuned': df[r['uid']]})
    if rows:
        with open(out / 'eval_init_vs_finetuned.csv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        try:
            from scipy.stats import wilcoxon
        except ImportError:
            wilcoxon = None
        for name in ('val', 'test'):
            sub = [r for r in rows if r['split'] == name]
            if not sub:
                continue
            log(f'--- {name}: Dice init -> fine-tuned (mean ± SD, n volumes) ---')
            groups = defaultdict(list)
            for r in sub:
                groups[f"{r['species']}/{r['dataset']}|{r['kind']}"].append(r)
            for k in sorted(groups):
                a = np.array([r['dice_init'] for r in groups[k]])
                b = np.array([r['dice_finetuned'] for r in groups[k]])
                log(f'    {k:<40} {a.mean():.4f}±{a.std():.4f} -> {b.mean():.4f}±{b.std():.4f} (n={len(a)})')
            a = np.array([r['dice_init'] for r in sub])
            b = np.array([r['dice_finetuned'] for r in sub])
            msg = f'    ALL {a.mean():.4f}±{a.std():.4f} -> {b.mean():.4f}±{b.std():.4f} (n={len(a)})'
            if wilcoxon is not None and len(a) >= 6 and np.any(a != b):
                st = wilcoxon(b, a)
                msg += f'  Wilcoxon signed-rank W={st.statistic:.1f}, p={st.pvalue:.2e}'
            log(msg)

    log(f'best checkpoint: {best_path} (epoch {best_epoch}, val macro-Dice {best:.4f})')
    log('to use it in EDNiX: copy it into ' + str(unet_dir / 'models') +
        " under a new name, then replace 'model-02-_cross_species-epoch' by that name in "
        'modalities/anat/skullstrip/Skullstrip_method.py and modalities/fMRI/Skullstrip_func.py')


if __name__ == '__main__':
    main()
