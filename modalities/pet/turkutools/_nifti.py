"""
Uncompressed NIfTI copies for the Turku PET Centre tools (tpcclib), which read .nii (not .nii.gz) and
take the frame times of a NIfTI image from a SIF with the same base name.

The PET engineer decompressed the input in place (gzip -d ... gzip -f): a failure in between left the
pipeline without its .nii.gz. Here the tools work on copies in a temporary folder next to the output
(inside the BIDS folder, which is bound in the Singularity container), and their outputs are written
back as .nii.gz with the geometry of the input image.
"""
import os
import shutil
import tempfile

import numpy as np
import nibabel as nib

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile


class TurkuWorkspace:
    """with TurkuWorkspace(out_dir) as ws: nii = ws.nii(img_gz, sif); ...; ws.back(nii_out, ref_gz, out_gz)"""

    def __init__(self, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        self.dir = tempfile.mkdtemp(prefix='tmp_turku_', dir=out_dir)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        shutil.rmtree(self.dir, ignore_errors=True)
        return False

    def path(self, name):
        return opj(self.dir, name)

    def nii(self, img_gz, sif=''):
        """Uncompressed copy of img_gz (and of its SIF, renamed to match)."""
        name = opb(img_gz)
        name = name[:-7] if name.endswith('.nii.gz') else name[:-4]
        out = self.path(name + '.nii')
        img = nib.load(img_gz)
        nib.save(nib.Nifti1Image(np.asarray(img.dataobj, dtype=np.float32), img.affine, img.header), out)
        if sif and opi(sif):
            shutil.copyfile(sif, self.path(name + '.sif'))
        return out

    @staticmethod
    def back(nii_out, ref_gz, out_gz):
        """Write the tool output nii_out as out_gz with the spatial geometry of ref_gz."""
        if not opi(nii_out):
            raise FileNotFoundError(f'Turku tool output not found: {nii_out}')
        ref = nib.load(ref_gz)
        out = nib.load(nii_out)
        data = np.asarray(out.dataobj, dtype=np.float32)
        if data.ndim > 3 and data.shape[3] == 1:
            data = data[..., 0]
        # a tool that wrote its output with another voxel order AND a valid header: bring it back to the
        # voxel order of the input (without a valid header the voxel order of the input is assumed)
        codes = (int(out.header['sform_code']), int(out.header['qform_code']))
        if max(codes) > 0 and not np.allclose(out.affine, ref.affine, atol=1e-2):
            o_out = nib.orientations.io_orientation(out.affine)
            o_ref = nib.orientations.io_orientation(ref.affine)
            if not np.array_equal(o_out, o_ref):
                data = nib.orientations.apply_orientation(data, nib.orientations.ornt_transform(o_out, o_ref))
        if data.shape[:3] != ref.shape[:3]:
            raise ValueError(f'{nii_out}: grid {data.shape[:3]} differs from the input {ref.shape[:3]}')
        hdr = ref.header.copy()
        hdr.set_data_dtype(np.float32)
        hdr.set_data_shape(data.shape)
        if data.ndim == 3 and len(hdr.get_zooms()) > 3:
            hdr.set_zooms(hdr.get_zooms()[:3])
        nib.save(nib.Nifti1Image(data, ref.affine, hdr), out_gz)
        return out_gz
