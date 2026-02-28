"""Analytical / numerical simulator for the two-variable ODE system.

Equations:
    dOmega/dt = -rho_tilde * Omega
    dChi/dt   = delta * Omega - phi * Chi

rho_tilde = rho_min + delta * (rho_max - rho_min)

This module exposes `simulate_ode(...)` for import and provides a CLI
for running simulations from the terminal and optionally saving/plotting
results.
"""

from typing import Tuple, Optional
import csv
import argparse
import sys

import numpy as np
from scipy.integrate import solve_ivp

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


def simulate_ode(
    delta: float,
    phi: float,
    rho_min: float,
    rho_max: float,
    t: np.ndarray,
    omega0: float = 1.0,
    chi0: float = 0.0,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    method: str = "DOP853",
) -> Tuple[np.ndarray, np.ndarray, float, Optional[float], Optional[float]]:
    """Numerically integrate the coupled ODEs and detect a peak of chi.

    Returns (Omega_t, Chi_t, rho_tilde, t_peak, chi_peak).
    """
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or t.size < 2:
        raise ValueError("t must be a 1D array with at least two time points")

    rho_tilde = rho_min + delta * (rho_max - rho_min)

    def rhs(_t, y):
        Omega, Chi = y
        dOmega = -rho_tilde * Omega
        dChi = delta * Omega - phi * Chi
        return [dOmega, dChi]

    def dchi_dt_event(_t, y):
        Omega, Chi = y
        return delta * Omega - phi * Chi

    # configure event: detect maxima where derivative crosses from + to -
    dchi_dt_event.terminal = False
    dchi_dt_event.direction = -1

    sol = solve_ivp(
        rhs,
        t_span=(t[0], t[-1]),
        y0=[omega0, chi0],
        t_eval=t,
        method=method,
        rtol=rtol,
        atol=atol,
        events=dchi_dt_event,
        dense_output=True,
    )

    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    Omega_t, Chi_t = sol.y

    # Peak detection: take first detected event (if any)
    if sol.t_events and len(sol.t_events) > 0 and sol.t_events[0].size > 0:
        t_peak = float(sol.t_events[0][0])
        _, Chi_peak = sol.sol(t_peak)
        chi_peak = float(Chi_peak)
    else:
        t_peak, chi_peak = None, None

    return Omega_t, Chi_t, float(rho_tilde), t_peak, chi_peak


def save_csv(path: str, t: np.ndarray, Omega: np.ndarray, Chi: np.ndarray) -> None:
    """Save time series to CSV with columns: t, Omega, Chi."""
    header = ["t", "Omega", "Chi"]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for ti, oi, ci in zip(t, Omega, Chi):
            writer.writerow([f"{ti:.10g}", f"{oi:.10g}", f"{ci:.10g}"])


def plot_results(t: np.ndarray, Omega: np.ndarray, Chi: np.ndarray, t_peak: Optional[float] = None, chi_peak: Optional[float] = None, save: Optional[str] = None) -> None:
    if plt is None:
        raise RuntimeError("matplotlib is not available for plotting")
    fig, ax = plt.subplots()
    ax.plot(t, Omega, label="Omega(t)")
    ax.plot(t, Chi, label="Chi(t)")
    if t_peak is not None and chi_peak is not None:
        ax.plot(t_peak, chi_peak, "ro", label="Chi peak")
    ax.set_xlabel("t")
    ax.set_ylabel("state")
    ax.legend()
    ax.grid(True)
    if save:
        fig.savefig(save, dpi=200)
    else:
        plt.show()


def main(argv=None):
    parser = argparse.ArgumentParser(description="Simulate Omega-Chi ODEs")
    parser.add_argument("--delta", type=float, default=0.5, help="coupling Δ in [0,1]")
    parser.add_argument("--phi", type=float, default=0.8, help="decay rate φ of χ")
    parser.add_argument("--rho-min", type=float, default=0.0)
    parser.add_argument("--rho-max", type=float, default=1.0)
    parser.add_argument("--t-max", type=float, default=10.0)
    parser.add_argument("--n-steps", type=int, default=500)
    parser.add_argument("--omega0", type=float, default=1.0)
    parser.add_argument("--chi0", type=float, default=0.0)
    parser.add_argument("--method", type=str, default="DOP853", help="integration method for solve_ivp")
    parser.add_argument("--out-csv", type=str, default=None, help="path to save CSV (optional)")
    parser.add_argument("--plot", action="store_true", help="show a plot (requires matplotlib)")
    parser.add_argument("--plot-save", type=str, default=None, help="save plot to file instead of showing")

    args = parser.parse_args(argv)

    t = np.linspace(0.0, args.t_max, args.n_steps)

    Omega, Chi, rho_tilde, t_peak, chi_peak = simulate_ode(
        delta=args.delta,
        phi=args.phi,
        rho_min=args.rho_min,
        rho_max=args.rho_max,
        t=t,
        omega0=args.omega0,
        chi0=args.chi0,
        method=args.method,
    )

    print(f"rho_tilde = {rho_tilde:.6g}")
    if t_peak is not None:
        print(f"chi peak at t = {t_peak:.6g}, chi = {chi_peak:.6g}")
    else:
        print("no chi peak detected in simulated interval")

    if args.out_csv:
        save_csv(args.out_csv, t, Omega, Chi)
        print(f"results saved to {args.out_csv}")

    if args.plot:
        try:
            plot_results(t, Omega, Chi, t_peak, chi_peak, save=args.plot_save)
            if args.plot_save:
                print(f"plot saved to {args.plot_save}")
        except Exception as exc:
            print(f"plot failed: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()