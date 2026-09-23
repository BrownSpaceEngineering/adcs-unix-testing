# ADCS Unix Testing

## Build 
```bash
git clone https://github.com/BrownSpaceEngineering/adcs-unix-testing.git
cd adcs-unix-testing
git submodule update --init --recursive
cmake -S . -B build
cmake --build build -j4
```

## Run the tests
```bash
./build/adcs-test
```
Every check prints `[PASS]` or `[FAIL]` (numeric checks also print the measured value and the
limit), and a summary at the end. A failing check never stops the run; the exit status is 1 if
anything failed. Tests live in `src/test.c`.

Some expected values come from `tools/reference_values.py`, a pure-Python (no numpy)
double-precision reference that re-derives them independently (e.g. the WMM field as
`-grad(V)` of the spherical-harmonic potential). Re-run it if you change a model:
```bash
python3 tools/reference_values.py
```

## Conventions
- Quaternions are Hamilton, scalar-first `[w, x, y, z]`, active (`quat_apply(q, v) = q v q^-1`).
  The filter estimates body->ECI.
- `body()` expects the magnetometer in body frame, **Tesla** (same units as the WMM model),
  the gyro in rad/s, and orbital elements as `[a (m), e, i, RAAN, argp, nu]` with angles in
  **degrees**. `orbital_to_eci*` take angles in radians.
- `iterate()` (the UKF) takes 0 (predict only), 1 (magnetometer only, eclipse) or 2 (magnetometer
  and sun) measurement vectors, which should be unit vectors.

## C vs MATLAB parity tests
The `MATLAB parity` suite in `src/test.c` replays inputs through the C code and compares against
outputs recorded from the team's MATLAB (`../adcs/MATLAB/Algorithms`) in
`include/matlab_reference.h`. The header is committed, so `adcs-test` never needs MATLAB. To
regenerate it after changing the MATLAB (runs MATLAB headless, then rebuilds and runs the tests):
```bash
tools/matlab/run_matlab_parity.sh            # or: tools/matlab/run_matlab_parity.sh /path/to/adcs
```
- Generation is deterministic (fixed seed), so the header only changes when the MATLAB does.
- Known MATLAB bugs are printed as `[NOTE]` lines, not failures. For those, the C code is compared
  with the MATLAB after one exact, named text fix applied to a temporary copy
  (`generate_matlab_reference.m`); nothing in `../adcs` is modified. If upstream MATLAB changes so
  a fix no longer applies, the header records it and the notes say so.
- `pointing_error.m` needs `quatmultiply`/`quatinv`/`rotm2quat` (Aerospace/Robotics toolboxes).
  If they're missing, the plain implementations in `tools/matlab/shims` are used.
- The UKF comparison uses the local functions of `simulink.m` directly (not `fcn`, which needs the
  Statistics Toolbox for `normrnd`).
