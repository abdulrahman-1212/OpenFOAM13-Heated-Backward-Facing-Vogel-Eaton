#!/usr/bin/env python3
"""
Compare CFD Nusselt number with Vogel & Eaton data (Ansys verification doc).

Reads bfs_results.csv (x/h, St, Cf, Twall) written by postprocessing.py and
converts St -> Nu_h = St * Re_h * Pr, with all properties taken from the
OpenFOAM case so the conversion is consistent with the simulation.

Usage:  python3 compareNu.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------- case properties
h = 0.0381          # step height [m]
U = 11.5            # reference velocity [m/s]
p = 1e5             # pressure [Pa]
Tin = 300.0         # inlet temperature [K]
W = 28.96           # molecular weight [kg/kmol]
mu = 1.85e-5        # dynamic viscosity [Pa s]  (physicalProperties)
Cp = 1005.0         # [J/kg/K]
Pr = 0.71

rho = p / (8314.462 / W * Tin)
k = mu * Cp / Pr
Re = rho * U * h / mu
print(f"rho = {rho:.4f} kg/m3, k = {k:.5f} W/m/K, Re_h = {Re:.0f}")

# ------------------------------------------------------------ CFD data
d = np.genfromtxt("bfs_results.csv", delimiter=",", skip_header=1)
x_cfd = d[:, 0] * h          # [m] from the step
St = d[:, 1]
Nu_cfd = St * Re * Pr        # equals q*h / (k*(Tw - Tin))

# ------------------------------------------------- experimental data
position = np.array([
    0.000, 0.034, 0.064, 0.095, 0.128,
    0.159, 0.190, 0.222, 0.253, 0.285,
    0.317, 0.350, 0.380, 0.411, 0.443,
    0.476, 0.507, 0.539, 0.571, 0.604
])
nusselt = np.array([
    0.0, 24.5, 29.2, 38.9, 48.7,
    55.5, 61.0, 62.8, 62.4, 60.7,
    58.5, 55.4, 52.4, 51.2, 49.9,
    48.6, 47.8, 47.4, 46.9, 46.5
])

# --------------------------------------------------------- comparison
Nu_at_exp = np.interp(position, x_cfd, Nu_cfd)
mask = position > 0
err = 100.0 * (Nu_at_exp[mask] - nusselt[mask]) / nusselt[mask]

print("\n  x [m]   x/h   Nu_exp   Nu_CFD   error %")
for xi, ne, nc in zip(position, nusselt, Nu_at_exp):
    e = "" if xi == 0 else f"{100*(nc-ne)/ne:7.1f}"
    print(f"  {xi:5.3f}  {xi/h:5.2f}  {ne:6.1f}  {nc:7.1f}  {e}")

print(f"\nMean abs error : {np.mean(np.abs(err)):.1f} %")
print(f"Max  abs error : {np.max(np.abs(err)):.1f} %")
print(f"Peak Nu  exp   : {nusselt.max():.1f} at x = {position[nusselt.argmax()]:.3f} m "
      f"(x/h = {position[nusselt.argmax()]/h:.2f})")
print(f"Peak Nu  CFD   : {Nu_cfd.max():.1f} at x = {x_cfd[Nu_cfd.argmax()]:.3f} m "
      f"(x/h = {x_cfd[Nu_cfd.argmax()]/h:.2f})")

# ---------------------------------------------------------------- plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x_cfd, Nu_cfd, "b-", lw=2, label="CFD (k-omega SST)")
ax.plot(position, nusselt, "ko", mfc="none", ms=7, label="Vogel & Eaton (exp.)")
ax.set_xlim(0, 0.62)
ax.set_xlabel("Distance from step, x [m]")
ax.set_ylabel("Nusselt number Nu")
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig("nu_comparison.png", dpi=150)
print("\nWrote nu_comparison.png")