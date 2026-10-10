#!/usr/bin/env python3
"""
Impose a turbulent boundary-layer inlet profile on the heated backward-facing step.

What it does
  * reads the inlet face centres from the C field written by
        foamPostProcess -latestTime -func writeCellCentres
  * builds profiles that match the Ansys VMFL013 inlet description:
        - 1/7 power-law velocity profile, boundary-layer thickness delta = 1.1 h
          on BOTH channel walls, uniform free stream Uinf set by Re_h = 28000
        - k   : free-stream turbulence + near-wall boundary-layer turbulence
        - omega: free stream from viscosity ratio, near wall from mixing length
  * replaces the 'inlet' block in 0/U, 0/k, 0/omega with nonuniform lists
    and sets the internal fields to the free-stream values
  * saves the original files as <name>.orig and plots the profiles

Usage (run from the case directory):
    python3 inletProfile.py                 # uses latest time that has a C file
    python3 inletProfile.py --time 5000
"""
import argparse
import os
import re
import shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--case", default=".")
ap.add_argument("--time", default=None, help="time dir holding the C file")
ap.add_argument("--field-dir", default="0", help="dir whose fields are edited")
ap.add_argument("--patch", default="inlet")
ap.add_argument("--h", type=float, default=0.0381, help="step height [m]")
ap.add_argument("--channel", type=float, default=4.0, help="inlet channel height / h")
ap.add_argument("--delta", type=float, default=1.1, help="BL thickness / h")
ap.add_argument("--Re", type=float, default=28000.0, help="Re_h on free-stream velocity")
ap.add_argument("--rho", type=float, default=1.161)
ap.add_argument("--mu", type=float, default=1.85e-5)
ap.add_argument("--I", type=float, default=0.02, help="free-stream turbulence intensity")
ap.add_argument("--nutRatio", type=float, default=10.0, help="free-stream nut/nu")
args = ap.parse_args()

Cmu, kappa = 0.09, 0.41
h = args.h
Hc = args.channel * h
delta = args.delta * h
nu = args.mu / args.rho
Uinf = args.Re * args.mu / (args.rho * h)


# ------------------------------------------------------------------ parsing
def find_block(text, name):
    bf = text.index("boundaryField")
    m = re.compile(r"\b" + re.escape(name) + r"\s*\{").search(text, bf)
    if not m:
        raise KeyError(f"patch '{name}' not found")
    i, depth = m.end(), 1
    while depth:
        c = text[i]
        depth += (c == "{") - (c == "}")
        i += 1
    return m.start(), i


def read_vector_patch(path, patch):
    text = open(path).read()
    s, e = find_block(text, patch)
    body = text[s:e]
    m = re.search(r"value\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", body)
    if not m:
        raise ValueError(f"no nonuniform vector list for {patch} in {path}")
    n = int(m.group(1))
    pat = r"\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)"
    rows = [g.groups() for g, _ in zip(re.finditer(pat, body[m.end():]), range(n))]
    return np.array(rows, dtype=float)


def latest_time_with_C(case):
    best = None
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isfile(os.path.join(case, d, "C")) and (best is None or t > best[0]):
            best = (t, d)
    if best is None:
        raise SystemExit("No time directory with a C file. Run:\n"
                         "  foamPostProcess -latestTime -func writeCellCentres")
    return best[1]


# ----------------------------------------------------------------- profiles
tdir = args.time or latest_time_with_C(args.case)
C = read_vector_patch(os.path.join(args.case, tdir, "C"), args.patch)
y = C[:, 1]
print(f"Read {len(y)} inlet faces from {tdir}/C   (y = {y.min():.2e} .. {y.max():.4f} m)")

d = np.minimum(y, Hc - y)                  # distance to nearest wall
eta = np.minimum(d / delta, 1.0)

U = Uinf * eta ** (1.0 / 7.0)

Re_d = Uinf * delta / nu
u_tau = Uinf * np.sqrt(0.0225 * Re_d ** -0.25)      # Prandtl 1/7-law friction

k_inf = 1.5 * (args.I * Uinf) ** 2
k = k_inf + (u_tau ** 2 / np.sqrt(Cmu)) * (1.0 - eta) ** 2

om_inf = k_inf / (args.nutRatio * nu)
l_mix = np.maximum(np.minimum(kappa * d, 0.09 * delta), 1e-8)
omega = np.maximum(om_inf, np.sqrt(k) / (Cmu ** 0.25 * l_mix))

# bulk velocity (info only)
o = np.argsort(y)
Ubulk = np.trapz(U[o], y[o]) / Hc
print(f"Uinf (free stream)  = {Uinf:.3f} m/s   (Re_h = {args.Re:.0f})")
print(f"Ubulk               = {Ubulk:.3f} m/s")
print(f"u_tau (estimate)    = {u_tau:.3f} m/s")
print(f"k_inf               = {k_inf:.4g} m2/s2,  omega_inf = {om_inf:.4g} 1/s")
print(f"max omega at inlet  = {omega.max():.3g} 1/s")


# ------------------------------------------------------------- file editing
def fmt_list(values, vector):
    if vector:
        lines = [f"({v:.6g} 0 0)" for v in values]
        typ = "vector"
    else:
        lines = [f"{v:.6g}" for v in values]
        typ = "scalar"
    return (f"List<{typ}> {len(values)}\n(\n" + "\n".join(lines) + "\n)")


def edit_field(name, values, vector, internal):
    path = os.path.join(args.case, args.field_dir, name)
    if not os.path.isfile(path):
        raise SystemExit("Missing " + path)
    if not os.path.isfile(path + ".orig"):
        shutil.copy(path, path + ".orig")
    text = open(path).read()

    s, e = find_block(text, args.patch)
    new = (f"{args.patch}\n    {{\n        type            fixedValue;\n"
           f"        value           nonuniform {fmt_list(values, vector)}\n;\n    }}")
    text = text[:s] + new + text[e:]

    text = re.sub(r"internalField\s+uniform\s+[^;]+;",
                  f"internalField   uniform {internal};", text, count=1)
    open(path, "w").write(text)
    print("Updated", path)


edit_field("U", U, True, f"({Uinf:.6g} 0 0)")
edit_field("k", k, False, f"{k_inf:.6g}")
edit_field("omega", omega, False, f"{om_inf:.6g}")

# ---------------------------------------------------------------------- plot
np.savetxt(os.path.join(args.case, "inlet_profile.csv"),
           np.column_stack([y[o] / h, U[o], k[o], omega[o]]),
           delimiter=",", header="y/h,U,k,omega", comments="")

fig, ax = plt.subplots(1, 3, figsize=(11, 4), sharey=True)
ax[0].plot(U[o], y[o] / h)
ax[0].set_xlabel("U [m/s]")
ax[0].set_ylabel("y/h")
ax[1].plot(k[o], y[o] / h)
ax[1].set_xlabel("k [m2/s2]")
ax[2].semilogx(omega[o], y[o] / h)
ax[2].set_xlabel("omega [1/s]")
for a in ax:
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(args.case, "inlet_profile.png"), dpi=130)
print("Wrote inlet_profile.png and inlet_profile.csv")