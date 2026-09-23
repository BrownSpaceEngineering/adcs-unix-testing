#!/usr/bin/env python3
"""
Independent double-precision reference values for src/test.c.

Pure Python (no numpy) so it runs anywhere:  python3 tools/reference_values.py

Each section re-derives the quantity a different way than the C code where
possible, so a shared bug can't hide:
  * WMM: B = -grad(V) of the spherical-harmonic potential by central finite
    differences in Cartesian ECEF, with Schmidt-normalised Legendre functions
    built from the explicit polynomial formula (the C code uses recursions,
    a geodetic round trip and a NED rotation instead).
  * Kepler: same algorithm as MATLAB/C but in double precision.
  * Sun vector / GMST: same formulas in double precision (float32 in C).

The coefficients are parsed from src/magnetosphere.c and cross-checked
against ../adcs/MATLAB/Algorithms/magnetosphere.m when that file exists.
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


# --------------------------------------------------------------------------
# WMM2025
# --------------------------------------------------------------------------
def load_c_coeffs():
    src = open(os.path.join(ROOT, "src", "magnetosphere.c")).read()
    rows = re.findall(
        r"\{\s*(\d+),\s*(\d+),\s*([-\d.]+)f,\s*([-\d.]+)f,\s*([-\d.]+)f,\s*([-\d.]+)f\s*\}", src)
    return [(int(n), int(m), float(g), float(h), float(gd), float(hd))
            for n, m, g, h, gd, hd in rows]


def load_matlab_coeffs():
    path = os.path.join(ROOT, "..", "adcs", "MATLAB", "Algorithms", "magnetosphere.m")
    if not os.path.exists(path):
        return None
    txt = open(path).read()
    block = txt[txt.index("data = [...") + len("data = [..."):txt.index("];", txt.index("data = [..."))]
    out = []
    for line in block.strip().splitlines():
        parts = line.split()
        if len(parts) == 6:
            n, m = int(parts[0]), int(parts[1])
            out.append((n, m) + tuple(float(p) for p in parts[2:]))
    return out


COEFFS = load_c_coeffs()
_ml = load_matlab_coeffs()
if _ml is not None:
    assert len(_ml) == len(COEFFS), "coefficient count mismatch vs MATLAB"
    for a, b in zip(COEFFS, _ml):
        assert a[:2] == b[:2] and all(abs(x - y) < 1e-9 for x, y in zip(a[2:], b[2:])), (a, b)

A_REF = 6371200.0


def jd2year(jd):
    return 2000.0 + (jd - 2451545.0) / 365.25


def gmst_wmm(jd):
    """GMST in radians, same polynomial as magnetosphere.m."""
    t = (jd - 2451545.0) / 36525.0
    sec = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * t + 0.093104 * t * t - 6.2e-6 * t ** 3
    return math.radians((sec / 240.0) % 360.0)


def gmst_ecef2eci(jd):
    """GMST in radians, same polynomial as eceftoeci.m."""
    t = (jd - 2451545.0) / 36525.0
    deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t * t - t ** 3 / 38710000.0
    return math.radians(deg % 360.0)


def legendre_poly_coeffs(n):
    """Coefficients c[k] of P_n(x) = sum c[k] x^k (explicit formula)."""
    c = [0.0] * (n + 1)
    for k in range(n // 2 + 1):
        c[n - 2 * k] = ((-1) ** k) * math.comb(n, k) * math.comb(2 * n - 2 * k, n) / 2.0 ** n
    return c


def schmidt_P(n, m, x):
    c = legendre_poly_coeffs(n)
    for _ in range(m):  # m-th derivative
        c = [k * c[k] for k in range(1, len(c))]
    val = sum(ck * x ** k for k, ck in enumerate(c))
    val *= (1.0 - x * x) ** (m / 2.0)
    if m > 0:
        val *= math.sqrt(2.0 * math.factorial(n - m) / math.factorial(n + m))
    return val


def potential(xyz, g, h):
    x, y, z = xyz
    r = math.sqrt(x * x + y * y + z * z)
    ct = z / r
    lon = math.atan2(y, x)
    v = 0.0
    for (n, m), gv in g.items():
        hv = h[(n, m)]
        v += (A_REF / r) ** (n + 1) * (gv * math.cos(m * lon) + hv * math.sin(m * lon)) * schmidt_P(n, m, ct)
    return A_REF * v


def wmm_eci(r_eci, jd):
    dy = jd2year(jd) - 2025.0
    g = {(n, m): gg + dy * gd for n, m, gg, hh, gd, hd in COEFFS}
    h = {(n, m): hh + dy * hd for n, m, gg, hh, gd, hd in COEFFS}
    th = gmst_wmm(jd)
    c, s = math.cos(th), math.sin(th)
    r_ecef = [c * r_eci[0] + s * r_eci[1], -s * r_eci[0] + c * r_eci[1], r_eci[2]]
    step = 10.0
    b_ecef = []
    for i in range(3):
        p = list(r_ecef)
        q = list(r_ecef)
        p[i] += step
        q[i] -= step
        b_ecef.append(-(potential(p, g, h) - potential(q, g, h)) / (2 * step))
    b_eci = [c * b_ecef[0] - s * b_ecef[1], s * b_ecef[0] + c * b_ecef[1], b_ecef[2]]
    return [b * 1e-9 for b in b_eci]  # nT -> T


# --------------------------------------------------------------------------
# Kepler / orbits
# --------------------------------------------------------------------------
MU = 3.986004418e14


def propagate(a, e, nu_deg, dt):
    f = math.radians(nu_deg)
    E = 2 * math.atan2(math.sqrt(1 - e) * math.sin(f / 2), math.sqrt(1 + e) * math.cos(f / 2))
    M = E - e * math.sin(E)
    M_new = M + math.sqrt(MU / a ** 3) * dt
    Ec = M_new
    for _ in range(100):
        d = (Ec - e * math.sin(Ec) - M_new) / (1 - e * math.cos(Ec))
        Ec -= d
        if abs(d) < 1e-15:
            break
    return math.degrees(2 * math.atan2(math.sqrt(1 + e) * math.sin(Ec / 2), math.sqrt(1 - e) * math.cos(Ec / 2)))


# --------------------------------------------------------------------------
# Sun
# --------------------------------------------------------------------------
def sun_vec(unix):
    jd = unix / 86400.0 + 2440587.5
    T = (jd - 2451545.0) / 36525.0
    M = (357.529 + 35999.050 * T) % 360
    L = (280.459 + 36000.770 * T) % 360
    C = ((1.914602 - 0.004817 * T - 0.000014 * T * T) * math.sin(math.radians(M))
         + (0.019993 - 0.000101 * T) * math.sin(math.radians(2 * M))
         + 0.000289 * math.sin(math.radians(3 * M)))
    lam = (L + C) % 360
    eps = 23 + 26 / 60 + 21.448 / 3600 - (46.8150 * T + 0.00059 * T * T - 0.001813 * T ** 3) / 3600
    return [math.cos(math.radians(lam)),
            math.cos(math.radians(eps)) * math.sin(math.radians(lam)),
            math.sin(math.radians(eps)) * math.sin(math.radians(lam))]


def fmt(v):
    return "{" + ", ".join("%.9gf" % x for x in v) + "}"


if __name__ == "__main__":
    print("// ---- WMM (Tesla, ECI) ----")
    cases = [
        ([6378137.0 + 500e3, 0.0, 0.0], 2461116, 0.25),
        ([0.0, -4.2e6, 5.4e6], 2461116, 0.25),
        ([3.1e6, 2.2e6, -5.6e6], 2461300, 0.7),
        ([-4.0e6, 1.0e6, 5.3e6], 2460900, 0.1),
        ([1.0e3, 2.0e3, 6.9e6], 2461000, 0.5),   # near north pole
    ]
    for r, jdi, jdf in cases:
        b = wmm_eci(r, jdi + jdf)
        print("{%s, %d, %.6ff, %s}," % (fmt(r), jdi, jdf, fmt(b)))

    print("// ---- GMST (rad) ----")
    for jd in (2451545.0, 2461116.25, 2461300.7):
        print("jd=%.6f  wmm=%.9f  ecef2eci=%.9f" % (jd, gmst_wmm(jd), gmst_ecef2eci(jd)))

    print("// ---- ecef2eci ----")
    unix = 1772814426
    jd = unix / 86400.0 + 2440587.5
    th = gmst_ecef2eci(jd)
    r = [6.3761e6, -0.1387e6, 0.0807e6]
    print("unix=%d r_eci=%s" % (unix, fmt([math.cos(th) * r[0] - math.sin(th) * r[1],
                                          math.sin(th) * r[0] + math.cos(th) * r[1], r[2]])))

    print("// ---- Kepler (new true anomaly, deg) ----")
    print("curtis: %.9f deg (%.6f rad)" % (propagate((9600e3 + 21000e3) / 2, 0.37255, 0.0, 10800.0),
                                           math.radians(propagate((9600e3 + 21000e3) / 2, 0.37255, 0.0, 10800.0)) % (2 * math.pi)))
    for a, e, nu, dt in ((7.0e6, 0.001, 10.0, 600.0), (1.0e7, 0.7, 1.0, 3600.0), (6.9e6, 0.05, 250.0, 1234.5)):
        print("a=%g e=%g nu=%g dt=%g -> %.9f" % (a, e, nu, dt, propagate(a, e, nu, dt)))

    print("// ---- Sun vector ----")
    for unix in (1710903960, 1718916660, 1767225600, 1790000000):
        print("unix=%d %s" % (unix, fmt(sun_vec(unix))))
