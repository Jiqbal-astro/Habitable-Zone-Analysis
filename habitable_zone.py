#!/usr/bin/env python3
"""
habitable_zone.py
=================

Habitable zone (HZ) boundaries for main-sequence stars following
Kopparapu et al. (2014), plus tools to place known exoplanets in their
star's habitable zone.

Features
--------
* Conservative HZ (runaway greenhouse -> maximum greenhouse) and
  optimistic HZ (recent Venus -> early Mars) limits
* Planet-mass dependent inner edge (0.1, 1 and 5 Earth masses)
* Insolation, equilibrium temperature and HZ distance (Kane & Gelino 2012)
  for any planet
* Orbital periods at the HZ edges from Kepler's third law
* Publication-style figures: the Teff-insolation HZ diagram, per-system
  orbit plots and a comparison of HZ extents

Author: Javed Iqbal (https://github.com/Jiqbal-astro)
"""

import argparse
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Kopparapu et al. (2014), Table 1
# S_eff = S_eff_sun + a*T + b*T^2 + c*T^3 + d*T^4,   T = Teff - 5780 K
# --------------------------------------------------------------------------
HZ_COEFFICIENTS = {
    "recent_venus":           (1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15),
    "runaway_greenhouse":     (1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15),
    "maximum_greenhouse":     (0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16),
    "early_mars":             (0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16),
    "runaway_greenhouse_5me": (1.188, 1.433e-4, 1.707e-8, -8.968e-12, -2.084e-15),
    "runaway_greenhouse_01me": (0.990, 1.209e-4, 1.404e-8, -7.418e-12, -1.713e-15),
}

# Inner-edge (runaway greenhouse) coefficient set for each planet mass
RUNAWAY_BY_MASS = {
    0.1: "runaway_greenhouse_01me",
    1.0: "runaway_greenhouse",
    5.0: "runaway_greenhouse_5me",
}

TEFF_VALID_RANGE = (2600.0, 7200.0)
DAYS_PER_YEAR = 365.25

DATA_DIR = Path(__file__).resolve().parent / "data"


# --------------------------------------------------------------------------
# Core physics
# --------------------------------------------------------------------------
def effective_flux(teff, limit):
    """Effective stellar flux S_eff (Earth = 1) at a given HZ limit."""
    if limit not in HZ_COEFFICIENTS:
        raise ValueError(f"Unknown limit '{limit}'. Choose from {list(HZ_COEFFICIENTS)}")
    s_sun, a, b, c, d = HZ_COEFFICIENTS[limit]
    t = np.asarray(teff, dtype=float) - 5780.0
    return s_sun + a * t + b * t**2 + c * t**3 + d * t**4


def check_teff(teff, name="star"):
    """Warn if Teff is outside the range where the Kopparapu fits are valid."""
    lo, hi = TEFF_VALID_RANGE
    if not lo <= teff <= hi:
        warnings.warn(
            f"{name}: Teff = {teff:.0f} K is outside the {lo:.0f}-{hi:.0f} K range of "
            "Kopparapu et al. (2014); HZ limits are extrapolated.",
            stacklevel=2,
        )
        return False
    return True


def hz_limits(teff, luminosity, planet_mass=1.0, name="star"):
    """Habitable zone boundaries in AU.

    Parameters
    ----------
    teff : float
        Stellar effective temperature (K).
    luminosity : float
        Stellar luminosity (solar units).
    planet_mass : {0.1, 1.0, 5.0}
        Planet mass in Earth masses; only affects the runaway greenhouse
        (conservative inner) edge.

    Returns
    -------
    dict with keys optimistic_inner, conservative_inner,
    conservative_outer, optimistic_outer (AU).
    """
    if planet_mass not in RUNAWAY_BY_MASS:
        raise ValueError(f"planet_mass must be one of {list(RUNAWAY_BY_MASS)}")
    check_teff(teff, name)

    def distance(limit):
        return float(np.sqrt(luminosity / effective_flux(teff, limit)))

    return {
        "optimistic_inner": distance("recent_venus"),
        "conservative_inner": distance(RUNAWAY_BY_MASS[planet_mass]),
        "conservative_outer": distance("maximum_greenhouse"),
        "optimistic_outer": distance("early_mars"),
    }


def insolation(luminosity, a):
    """Stellar flux received by a planet, in units of Earth's insolation."""
    return luminosity / np.asarray(a, dtype=float) ** 2


def equilibrium_temperature(luminosity, a, albedo=0.3):
    """Planetary equilibrium temperature (K), assuming full heat redistribution.

    Earth (L = 1, a = 1 AU, albedo 0.3) gives ~255 K.
    """
    return 278.6 * (1.0 - albedo) ** 0.25 * luminosity ** 0.25 / np.sqrt(a)


def hz_distance(a, inner, outer):
    """Habitable Zone Distance (Kane & Gelino 2012).

    -1 at the inner edge, 0 at the centre, +1 at the outer edge.
    Values between -1 and +1 lie inside the HZ.
    """
    return (2.0 * a - outer - inner) / (outer - inner)


def orbital_period_days(a, stellar_mass):
    """Orbital period (days) from Kepler's third law, neglecting planet mass."""
    return DAYS_PER_YEAR * np.sqrt(np.asarray(a, dtype=float) ** 3 / stellar_mass)


def classify(a, limits):
    """Where a planet sits relative to the HZ."""
    if limits["conservative_inner"] <= a <= limits["conservative_outer"]:
        return "conservative HZ"
    if limits["optimistic_inner"] <= a <= limits["optimistic_outer"]:
        return "optimistic HZ"
    return "too hot" if a < limits["optimistic_inner"] else "too cold"


# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------
def load_catalogue(stars_csv=DATA_DIR / "stars.csv", planets_csv=DATA_DIR / "planets.csv"):
    stars = pd.read_csv(stars_csv)
    planets = pd.read_csv(planets_csv) if Path(planets_csv).exists() else pd.DataFrame(
        columns=["planet", "star", "semi_major_axis_au"])
    return stars, planets


def star_table(stars, planet_mass=1.0):
    """HZ boundaries and edge orbital periods for every star."""
    rows = []
    for _, s in stars.iterrows():
        lim = hz_limits(s.teff_k, s.luminosity_lsun, planet_mass, name=s.star)
        rows.append({
            "star": s.star,
            "teff_k": s.teff_k,
            "luminosity_lsun": s.luminosity_lsun,
            "mass_msun": s.mass_msun,
            **{f"{k}_au": v for k, v in lim.items()},
            "hz_width_au": lim["conservative_outer"] - lim["conservative_inner"],
            "period_inner_days": orbital_period_days(lim["conservative_inner"], s.mass_msun),
            "period_outer_days": orbital_period_days(lim["conservative_outer"], s.mass_msun),
            "teff_in_valid_range": TEFF_VALID_RANGE[0] <= s.teff_k <= TEFF_VALID_RANGE[1],
        })
    return pd.DataFrame(rows)


def planet_table(stars, planets, planet_mass=1.0, albedo=0.3):
    """Insolation, temperature, HZ distance and classification for each planet."""
    star_info = stars.set_index("star")
    rows = []
    for _, p in planets.iterrows():
        if p.star not in star_info.index:
            warnings.warn(f"Host star '{p.star}' of {p.planet} not found; skipping.")
            continue
        s = star_info.loc[p.star]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lim = hz_limits(s.teff_k, s.luminosity_lsun, planet_mass)
        a = p.semi_major_axis_au
        rows.append({
            "planet": p.planet,
            "star": p.star,
            "semi_major_axis_au": a,
            "period_days": orbital_period_days(a, s.mass_msun),
            "insolation_searth": insolation(s.luminosity_lsun, a),
            "t_eq_k": equilibrium_temperature(s.luminosity_lsun, a, albedo),
            "hz_distance": hz_distance(a, lim["conservative_inner"], lim["conservative_outer"]),
            "location": classify(a, lim),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
LOCATION_COLOURS = {
    "conservative HZ": "#1a9850",
    "optimistic HZ": "#91cf60",
    "too hot": "#d73027",
    "too cold": "#4575b4",
}


def plot_hz_diagram(stars, planets_df, out_path, planet_mass=1.0):
    """Classic HZ diagram: stellar Teff vs insolation, with planets overplotted."""
    teff = np.linspace(*TEFF_VALID_RANGE, 400)
    s_rv = effective_flux(teff, "recent_venus")
    s_rg = effective_flux(teff, RUNAWAY_BY_MASS[planet_mass])
    s_mg = effective_flux(teff, "maximum_greenhouse")
    s_em = effective_flux(teff, "early_mars")

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.fill_betweenx(teff, s_em, s_rv, color="#91cf60", alpha=0.35, label="Optimistic HZ")
    ax.fill_betweenx(teff, s_mg, s_rg, color="#1a9850", alpha=0.45, label="Conservative HZ")
    for s, ls, lab in ((s_rv, ":", "Recent Venus"), (s_rg, "-", "Runaway greenhouse"),
                       (s_mg, "-", "Maximum greenhouse"), (s_em, ":", "Early Mars")):
        ax.plot(s, teff, color="k", ls=ls, lw=0.8)

    teff_by_star = stars.set_index("star")["teff_k"]
    for _, p in planets_df.iterrows():
        t = teff_by_star[p.star]
        ax.plot(p.insolation_searth, t, "o", ms=8, mec="k", mew=0.6,
                color=LOCATION_COLOURS[p.location], zorder=5)
        ax.annotate(p.planet, (p.insolation_searth, t), xytext=(6, 4),
                    textcoords="offset points", fontsize=8)

    ax.set_xscale("log")
    ax.set_xlim(4, 0.1)  # hot planets on the left, as is conventional
    ticks = [4, 2, 1, 0.5, 0.25, 0.1]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(2400, 6200)
    ax.axhspan(2400, TEFF_VALID_RANGE[0], color="0.85", alpha=0.6, lw=0)
    ax.text(3.6, 2480, "Teff < 2600 K: outside fit range", fontsize=8, color="0.35")
    ax.set_xlabel("Insolation (Earth = 1)")
    ax.set_ylabel("Stellar effective temperature (K)")
    ax.set_title("Habitable zone diagram (Kopparapu et al. 2014)")

    handles, labels = ax.get_legend_handles_labels()
    for loc, col in LOCATION_COLOURS.items():
        handles.append(plt.Line2D([], [], marker="o", ls="", color=col, mec="k", mew=0.6))
        labels.append(f"Planet: {loc}")
    ax.legend(handles, labels, loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_systems(stars_df, planets_df, out_path):
    """Top-down view of each system, each scaled to its own habitable zone."""
    n = len(stars_df)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).ravel()
    theta = np.linspace(0, 2 * np.pi, 300)

    for ax, (_, s) in zip(axes, stars_df.iterrows()):
        def ring(r_in, r_out, colour, alpha):
            ax.fill(np.r_[r_out * np.cos(theta), r_in * np.cos(theta[::-1])],
                    np.r_[r_out * np.sin(theta), r_in * np.sin(theta[::-1])],
                    color=colour, alpha=alpha, lw=0)

        ring(s.optimistic_inner_au, s.optimistic_outer_au, "#91cf60", 0.35)
        ring(s.conservative_inner_au, s.conservative_outer_au, "#1a9850", 0.45)

        sys_planets = planets_df[planets_df.star == s.star]
        extent = max([s.optimistic_outer_au] + list(sys_planets.semi_major_axis_au)) * 1.15
        for i, (_, p) in enumerate(sys_planets.iterrows()):
            a = p.semi_major_axis_au
            ax.plot(a * np.cos(theta), a * np.sin(theta), color="0.4", lw=0.6, ls="--")
            ang = np.deg2rad(35 + 55 * i)
            ax.plot(a * np.cos(ang), a * np.sin(ang), "o", ms=6, mec="k", mew=0.5,
                    color=LOCATION_COLOURS[p.location])
            ax.annotate(p.planet.replace(s.star, "").strip() or p.planet,
                        (a * np.cos(ang), a * np.sin(ang)), xytext=(4, 4),
                        textcoords="offset points", fontsize=8)
        ax.plot(0, 0, "*", ms=14, color="#fdb813", mec="#b8860b")
        ax.set_xlim(-extent, extent)
        ax.set_ylim(-extent, extent)
        ax.set_aspect("equal")
        ax.set_title(f"{s.star}\n{s.teff_k:.0f} K, HZ {s.conservative_inner_au:.3g}–"
                     f"{s.conservative_outer_au:.3g} AU", fontsize=10)
        ax.set_xlabel("AU", fontsize=8)
        ax.tick_params(labelsize=7)

    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Planetary systems and their habitable zones (each panel to its own scale)",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_hz_comparison(stars_df, planets_df, out_path):
    """Horizontal bars of HZ extent per star on a common log distance axis."""
    df = stars_df.sort_values("teff_k").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10, 0.6 * len(df) + 1.5))
    for i, s in df.iterrows():
        ax.barh(i, s.optimistic_outer_au - s.optimistic_inner_au, left=s.optimistic_inner_au,
                height=0.6, color="#91cf60", alpha=0.5)
        ax.barh(i, s.conservative_outer_au - s.conservative_inner_au,
                left=s.conservative_inner_au, height=0.6, color="#1a9850", alpha=0.8)
        sys_planets = planets_df[planets_df.star == s.star]
        ax.scatter(sys_planets.semi_major_axis_au, np.full(len(sys_planets), i), s=30,
                   c=[LOCATION_COLOURS[l] for l in sys_planets.location],
                   edgecolors="k", linewidths=0.5, zorder=5)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{s.star} ({s.teff_k:.0f} K)" for _, s in df.iterrows()])
    ax.set_xscale("log")
    lo = min(df.optimistic_inner_au.min(), planets_df.semi_major_axis_au.min()) * 0.8
    hi = max(df.optimistic_outer_au.max(), planets_df.semi_major_axis_au.max()) * 1.25
    ax.set_xlim(lo, hi)
    ticks = [t for t in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5) if lo <= t <= hi]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("Distance from star (AU)")
    ax.set_title("Habitable zone extent: conservative (dark) and optimistic (light)")
    ax.grid(alpha=0.3, axis="x", which="both")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------
def analyse_custom_star(args):
    lim = hz_limits(args.teff, args.lum, args.planet_mass, name=args.name)
    print(f"\nHabitable zone of {args.name} "
          f"(Teff = {args.teff:.0f} K, L = {args.lum:g} L_sun, M = {args.mass:g} M_sun)\n")
    for key in ("optimistic_inner", "conservative_inner", "conservative_outer", "optimistic_outer"):
        a = lim[key]
        print(f"  {key.replace('_', ' '):<20s} {a:8.4f} AU   "
              f"P = {orbital_period_days(a, args.mass):8.2f} d   "
              f"T_eq = {equilibrium_temperature(args.lum, a, args.albedo):5.0f} K")
    for a in args.planet_a or []:
        print(f"\n  planet at {a} AU: S = {insolation(args.lum, a):.2f} S_earth, "
              f"T_eq = {equilibrium_temperature(args.lum, a, args.albedo):.0f} K, "
              f"HZD = {hz_distance(a, lim['conservative_inner'], lim['conservative_outer']):+.2f}, "
              f"{classify(a, lim)}")


def run_catalogue(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    stars, planets = load_catalogue(args.stars, args.planets)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        stars_df = star_table(stars, args.planet_mass)
    for w in caught:
        print(f"Warning: {w.message}")
    planets_df = planet_table(stars, planets, args.planet_mass, args.albedo)

    stars_df.to_csv(out / "hz_stars.csv", index=False, float_format="%.5g")
    planets_df.to_csv(out / "hz_planets.csv", index=False, float_format="%.4g")

    pd.set_option("display.width", 140)
    print("\nHabitable zone limits (AU)")
    print(stars_df[["star", "teff_k", "optimistic_inner_au", "conservative_inner_au",
                    "conservative_outer_au", "optimistic_outer_au"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\nPlanets")
    print(planets_df[["planet", "semi_major_axis_au", "insolation_searth", "t_eq_k",
                      "hz_distance", "location"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    plot_hz_diagram(stars, planets_df, out / "hz_diagram.png", args.planet_mass)
    plot_systems(stars_df, planets_df, out / "hz_systems.png")
    plot_hz_comparison(stars_df, planets_df, out / "hz_comparison.png")
    print(f"\nTables and figures saved to {out}/")


def parse_args():
    p = argparse.ArgumentParser(
        description="Habitable zones (Kopparapu et al. 2014) for a catalogue of stars, "
                    "or for a single star given on the command line.")
    p.add_argument("-o", "--output", default="output", help="output directory (default: output)")
    p.add_argument("--stars", default=DATA_DIR / "stars.csv", help="star catalogue CSV")
    p.add_argument("--planets", default=DATA_DIR / "planets.csv", help="planet catalogue CSV")
    p.add_argument("--planet-mass", type=float, default=1.0, choices=[0.1, 1.0, 5.0],
                   help="planet mass in Earth masses for the inner edge (default: 1)")
    p.add_argument("--albedo", type=float, default=0.3, help="Bond albedo for T_eq (default: 0.3)")

    c = p.add_argument_group("single star (skips the catalogue)")
    c.add_argument("--teff", type=float, help="stellar effective temperature (K)")
    c.add_argument("--lum", type=float, help="stellar luminosity (L_sun)")
    c.add_argument("--mass", type=float, default=1.0, help="stellar mass (M_sun), for periods")
    c.add_argument("--name", default="Custom star", help="star name")
    c.add_argument("--planet-a", type=float, nargs="+", help="planet semi-major axes (AU) to assess")
    return p.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    if arguments.teff is not None and arguments.lum is not None:
        analyse_custom_star(arguments)
    else:
        run_catalogue(arguments)
