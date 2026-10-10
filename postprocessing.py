#!/usr/bin/env python3
"""
Post-processing for the heated backward-facing step (OpenFOAM 13, fluid module).

Needs, in the latest time directory (created by runPost.sh):
    T                  (heatedWall values)
    wallShearStress    (heatedWall values)
    C                  (cell centres; boundary values are face centres)

Outputs:
    bfs_results.csv    x/h, Stanton number, Cf, wall temperature
    bfs_results.png    St and Cf plots
    printed summary    reattachment length, St peak, mean y+ hint

Usage:
    python3 bfsPost.py [--case .] [--exp vogel_eaton.csv]

--exp is optional: a CSV with two columns "x/h, St" (header allowed) to
overlay your own digitised experimental data.
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- arguments
ap = argparse.ArgumentParser()
ap.add_argument("--case", default=".")
ap.add_argument("--patch", default="heatedWall")
ap.add_argument("--h", type=float, default=0.0381, help="step height [m]")
ap.add_argument("--q", type=float, default=270.0, help="wall heat flux [W/m2]")
ap.add_argument("--Tin", type=float, default=300.0, help="inlet temperature [K]")
ap.add_argument("--U", type=float, default=11.5, help="reference velocity [m/s]")
ap.add_argument("--Cp", type=float, default=1005.0)
ap.add_argument("--p", type=float, default=1e5, help="reference pressure [Pa]")
ap.add_argument("--W", type=float, default=28.96, help="molecular weight [kg/kmol]")
ap.add_argument("--exp", default=None, help="CSV with x/h, St (optional)")
args = ap.parse_args()

# ------------------------------------------------------------------ helpers
def latest_time(case):
    times = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(case, d)):
            times.append((t, d))
    if not times:
        sys.exit("No time directories found in " + case)
    return max(times)[1]


def block(text, name):
    """Return the {...} body that follows 'name' inside boundaryField."""
    bf = text.index("boundaryField")
    m = re.compile(r"\b" + re.escape(name) + r"\s*\{").search(text, bf)
    if not m:
        raise KeyError(f"patch '{name}' not found")
    i, depth = m.end(), 1
    while depth:
        c = text[i]
        depth += (c == "{") - (c == "}")
        i += 1
    return text[m.end():i - 1]


def read_patch(path, patch):
    """Read the 'value' list of one patch. Returns (array, dimensions list)."""
    text = open(path).read()
    dims = re.search(r"dimensions\s*\[([^\]]+)\]", text)
    dims = [float(v) for v in dims.group(1).split()] if dims else None
    body = block(text, patch)
    m = re.search(r"value\s+uniform\s+([^;]+);", body)
    if m:
        s = m.group(1).replace("(", " ").replace(")", " ").split()
        return np.array([float(v) for v in s]), dims
    m = re.search(r"value\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", body)
    if not m:
        raise ValueError(f"cannot parse value of patch {patch} in {path}")
    typ, n, pos = m.group(1), int(m.group(2)), m.end()
    rest = body[pos:]
    if typ == "scalar":
        vals = rest[:rest.index(")")].split()
        return np.array(vals[:n], dtype=float), dims
    if typ == "vector":
        pat = r"\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)"
        vecs = [m2.groups() for m2, _ in zip(re.finditer(pat, rest), range(n))]
        return np.array(vecs, dtype=float), dims
    raise ValueError("unsupported type " + typ)


# -------------------------------------------------------------------- read
t_dir = os.path.join(args.case, latest_time(args.case))
print("Reading time directory:", t_dir)

C, _ = read_patch(os.path.join(t_dir, "C"), args.patch)
T, _ = read_patch(os.path.join(t_dir, "T"), args.patch)
tau, tdims = read_patch(os.path.join(t_dir, "wallShearStress"), args.patch)

x = C[:, 0] / args.h
order = np.argsort(x)
x, T, tau_x = x[order], T[order], tau[order, 0]

R = 8314.462 / args.W
rho = args.p / (R * args.Tin)
print(f"Reference density (ideal gas, Tin={args.Tin} K): {rho:.4f} kg/m3")

# wallShearStress is in Pa for compressible solvers; convert if kinematic
if tdims is not None and abs(tdims[0]) < 1e-9:
    tau_x = tau_x * rho

# Orient sign so that attached flow (downstream end) is positive
if np.mean(tau_x[-max(5, len(tau_x) // 5):]) < 0:
    tau_x = -tau_x

Cf = tau_x / (0.5 * rho * args.U**2)
St = args.q / (rho * args.Cp * args.U * (T - args.Tin))

# ------------------------------------------------------------ reattachment
# Corner eddy gives tau_x>0 near the step, then the primary bubble gives
# tau_x<0. Reattachment is the last negative->positive crossing.
xr = None
for i in range(len(x) - 1):
    if x[i] > 0.0 and tau_x[i] < 0.0 <= tau_x[i + 1]:
        xr = x[i] + (x[i + 1] - x[i]) * (-tau_x[i]) / (tau_x[i + 1] - tau_x[i])

i_pk = int(np.argmax(St))
print("-" * 52)
if xr is None:
    print("Reattachment point: not found (no tau_x sign change; "
          "run longer or check convergence)")
else:
    print(f"Reattachment length  x_r/h = {xr:.2f}   (experiment ~ 6.7)")
print(f"St peak              {St[i_pk]:.5f} at x/h = {x[i_pk]:.2f}")
print(f"Wall T range         {T.min():.2f} .. {T.max():.2f} K")
print("-" * 52)

# ------------------------------------------------------------------ output
out_csv = os.path.join(args.case, "bfs_results.csv")
np.savetxt(out_csv, np.column_stack([x, St, Cf, T]), delimiter=",",
           header="x/h,St,Cf,Twall", comments="")
print("Wrote", out_csv)

fig, ax = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
ax[0].plot(x, St, "b-", label="CFD (k-omega SST)")
if args.exp:
    try:
        e = np.genfromtxt(args.exp, delimiter=",", skip_header=0)
        e = e[~np.isnan(e).any(axis=1)]
        ax[0].plot(e[:, 0], e[:, 1], "ko", mfc="none", label="Experiment")
    except Exception as err:
        print("Could not read experimental file:", err)
ax[0].set_ylabel("Stanton number St")
ax[0].legend()
ax[0].grid(alpha=0.3)

ax[1].plot(x, Cf, "r-")
ax[1].axhline(0, color="k", lw=0.8)
if xr is not None:
    ax[1].axvline(xr, color="g", ls="--", label=f"reattachment x/h = {xr:.2f}")
    ax[1].legend()
ax[1].set_xlabel("x/h")
ax[1].set_ylabel("Skin friction Cf")
ax[1].grid(alpha=0.3)

fig.tight_layout()
out_png = os.path.join(args.case, "bfs_results.png")
fig.savefig(out_png, dpi=150)
print("Wrote", out_png)