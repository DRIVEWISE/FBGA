"""
Tests for the fbga_py bindings to the fbga3d module (FbgaIndy/FbgaMoto), mirroring
test/fbga3d/test_fbga_solver.cc's flat-track scenario.
"""

import math

import pytest

import fbga_py as fb


def make_flat_trajectory():
    """Flat, non-banked track (mu = phi = 0 everywhere) with a simple curvature profile
    via theta_prime, mirroring test/fbga3d/test_fbga_solver.cc's make_flat_trajectory()."""
    TOA = fb.TrajectoryOffsetAndAnglesContainer()
    TOA.reference = fb.RoadAnglesAndDerivativesContainer(
        mu=[0.0] * 5,
        phi=[0.0] * 5,
        theta=[0.0, 0.05, 0.55, 0.8, 0.8],
        mu_prime=[0.0] * 5,
        phi_prime=[0.0] * 5,
        theta_prime=[0.001, 0.01, 0.005, -0.01, 0.0],
        abscissa=[0.0, 50.0, 100.0, 150.0, 200.0],
    )
    TOA.offset = fb.TrajectoryOffsetContainer(n=[0.0] * 5, chi=[0.0] * 5)
    TOA.adherence = fb.AdherenceContainer(alpha=[1.0] * 5)
    return TOA


def test_fbga_moto_compute_returns_finite_total_time():
    solver = fb.FbgaMoto()
    TOA = make_flat_trajectory()

    total_time = solver.compute(TOA, 20.0)

    assert math.isfinite(total_time)
    assert total_time > 0.0

    sol = solver.get_solution()
    assert len(sol.V0) == len(TOA.reference.abscissa) - 1
    assert all(math.isfinite(v) for v in sol.V0)
    assert all(math.isfinite(v) for v in sol.V1)


def test_fbga_moto_accessors_are_finite_along_the_solved_profile():
    solver = fb.FbgaMoto()
    TOA = make_flat_trajectory()
    solver.compute(TOA, 20.0)

    s = 1.0
    while s < 199.0:
        assert math.isfinite(solver.eval_V(s))
        assert math.isfinite(solver.eval_A_tilde_x(s))
        assert math.isfinite(solver.eval_Vmax(s))
        s += 20.0


def test_fbga_moto_segment_type_is_a_segment_type_enum():
    solver = fb.FbgaMoto()
    TOA = make_flat_trajectory()
    solver.compute(TOA, 20.0)

    assert isinstance(solver.eval_segment_type(10.0), fb.SegmentType)


def test_gggv_moto_model_is_reachable_and_settable_via_the_solver():
    solver = fb.FbgaMoto()
    # model() returns a mutable reference into the solver's own GggvMoto instance.
    assert solver.model().a_x_eng(20.0) > 0.0


def test_gggv_moto_a_x_push_shrinks_with_alpha():
    gggv = fb.GggvMoto()
    a_x_push_full = gggv.a_x_push(0.0, 20.0, fb.GRAVITY, 1.0)
    a_x_push_half = gggv.a_x_push(0.0, 20.0, fb.GRAVITY, 0.5)
    assert a_x_push_half <= a_x_push_full


def test_fbga_indy_compute_returns_finite_total_time():
    # GggvIndy.setup_std() (invoked by FbgaIndy's default constructor) reads
    # "./data/INDY/*.npy" relative to the CWD -- run pytest from the FBGA repo root.
    solver = fb.FbgaIndy()
    TOA = make_flat_trajectory()

    total_time = solver.compute(TOA, 20.0)

    assert math.isfinite(total_time)
    assert total_time > 0.0
