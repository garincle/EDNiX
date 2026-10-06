# Tests

## Unit tests

```bash
pip install pytest pytest-xdist
pytest -n auto
```

`tests/unit` (about 230 tests, a few minutes) checks, on synthetic data:

- that every Python file compiles and every pipeline module imports;
- geometry: orientation mappings (always proper rotations), grid comparisons, header checks;
- registration parameters of each step;
- functional preprocessing: normalisation, atlas erosion, one-shot resampling of motion, distortion correction (sign detection, field direction, order of the transforms, choice of the field per run);
- statistics: random-intercept maximum likelihood against statsmodels, type-I error of the longitudinal age test, cluster extent, component matching, phylogenetic covariance;
- launchers: every example launcher only uses parameters known to its pipeline, and the interface rewrites launchers without changing them;
- interface: label tables, quality-control files, run progress.

`pytest -m "not slow"` runs the fast subset in a few seconds. The tests run on every push (`.github/workflows/tests.yml`).

## End-to-end smoke test

On a machine with the data and the containers:

```bash
python tests/e2e/run_smoke.py                  # every example launcher whose data is on this machine
python tests/e2e/run_smoke.py --only anat
python tests/e2e/run_smoke.py --dry-run        # prepare the reduced data and launchers only
```

For each dataset it copies one subject (functional runs shortened to 40 volumes), writes a copy of the launcher pointing to that copy, skips the interactive, MELODIC and FreeSurfer / surface steps, and runs it with `EDNIX_STRICT=1`, anatomical launchers first. The report is written to `smoke_report.md` in the work folder (`/tmp/ednix_smoke` by default); the exit code is 1 if a run failed.

`.github/workflows/smoke.yml` runs it on demand on a self-hosted GitHub runner with the label `ednix`, installed on the server that holds the data and the containers.

## Parameter reference

When you change the signature or the docstring of a pipeline entry point, regenerate the reference pages:

```bash
python docs/generate_parameter_reference.py
```
