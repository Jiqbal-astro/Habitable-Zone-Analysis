"""Checks against published values. Run with:  pytest"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import habitable_zone as hz  # noqa: E402


def test_solar_system_conservative_hz():
    # Kopparapu et al. (2014): 0.95 - 1.67 AU for a 1 Earth-mass planet
    lim = hz.hz_limits(5780, 1.0)
    assert lim["conservative_inner"] == pytest.approx(0.95, abs=0.005)
    assert lim["conservative_outer"] == pytest.approx(1.676, abs=0.005)


def test_solar_system_optimistic_hz():
    lim = hz.hz_limits(5780, 1.0)
    assert lim["optimistic_inner"] == pytest.approx(0.75, abs=0.005)
    assert lim["optimistic_outer"] == pytest.approx(1.77, abs=0.005)


def test_bigger_planets_have_closer_inner_edge():
    small = hz.hz_limits(5780, 1.0, planet_mass=0.1)["conservative_inner"]
    earth = hz.hz_limits(5780, 1.0, planet_mass=1.0)["conservative_inner"]
    big = hz.hz_limits(5780, 1.0, planet_mass=5.0)["conservative_inner"]
    assert big < earth < small


def test_limits_are_finite_and_ordered_across_valid_range():
    for teff in np.linspace(2600, 7200, 50):
        lim = hz.hz_limits(teff, 1.0)
        values = [lim["optimistic_inner"], lim["conservative_inner"],
                  lim["conservative_outer"], lim["optimistic_outer"]]
        assert np.all(np.isfinite(values))
        assert values == sorted(values)


def test_warning_outside_valid_teff():
    with pytest.warns(UserWarning):
        hz.hz_limits(2500, 0.0005)


def test_earth_equilibrium_temperature():
    assert hz.equilibrium_temperature(1.0, 1.0, albedo=0.3) == pytest.approx(255, abs=1)


def test_earth_orbital_period():
    assert hz.orbital_period_days(1.0, 1.0) == pytest.approx(365.25)


def test_hz_distance_edges():
    assert hz.hz_distance(1.0, 1.0, 2.0) == pytest.approx(-1)
    assert hz.hz_distance(1.5, 1.0, 2.0) == pytest.approx(0)
    assert hz.hz_distance(2.0, 1.0, 2.0) == pytest.approx(1)


def test_classification_of_solar_system():
    lim = hz.hz_limits(5780, 1.0)
    assert hz.classify(0.723, lim) == "too hot"            # Venus
    assert hz.classify(1.0, lim) == "conservative HZ"      # Earth
    assert hz.classify(5.2, lim) == "too cold"             # Jupiter
