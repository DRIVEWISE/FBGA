"""
Tests for the fbga_py bindings, mirroring examples/fbga/basic.cc
(friction-circle G-G bounds driving the Fb2d forward-backward solver).
"""

import math

import pytest

import fbga_py as fb


MU_X = 1.3   # Longitudinal friction coefficient
MU_Y = 1.4   # Lateral friction coefficient
G = fb.GRAVITY


def gg_upper(ay: float, v: float) -> float:
    """Friction-circle upper bound on ax."""
    ay_norm = ay / G
    return G * math.sqrt(max(0.0, MU_Y * MU_Y - ay_norm * ay_norm)) * MU_X / MU_Y


def gg_lower(ay: float, v: float) -> float:
    """Friction-circle lower bound on ax."""
    return -gg_upper(ay, v)


def gg_range_min(v: float) -> float:
    return -MU_Y * G


def gg_range_max(v: float) -> float:
    return MU_Y * G


@pytest.fixture
def gg_range():
    r = fb.GgRangeMaxMin()
    r.min = gg_range_min
    r.max = gg_range_max
    return r


@pytest.fixture
def solver(gg_range):
    return fb.Fb2d(gg_upper, gg_lower, gg_range)


@pytest.fixture
def track():
    SS_vec = [0.0, 50.0, 100.0, 150.0, 200.0]
    KK_vec = [0.001, 0.01, 0.005, -0.01, 0.0]
    n_eval = math.ceil(SS_vec[-1])
    SS_eval = [float(i) for i in range(n_eval)]
    KK_eval = []
    for s in SS_eval:
        for j in range(len(SS_vec) - 1):
            if SS_vec[j] <= s <= SS_vec[j + 1]:
                t = (s - SS_vec[j]) / (SS_vec[j + 1] - SS_vec[j])
                KK_eval.append(KK_vec[j] + t * (KK_vec[j + 1] - KK_vec[j]))
                break
        else:
            KK_eval.append(KK_vec[-1])
    return SS_eval, KK_eval


def test_constants_exposed():
    assert fb.GRAVITY == pytest.approx(9.81)
    assert fb.PI == pytest.approx(math.pi)
    assert fb.DEG2RAD == pytest.approx(math.pi / 180.0)
    assert fb.RAD2DEG == pytest.approx(180.0 / math.pi)


def test_gg_range_roundtrip(gg_range):
    assert gg_range.min(10.0) == pytest.approx(-MU_Y * G)
    assert gg_range.max(10.0) == pytest.approx(MU_Y * G)


def test_fb2d_construction(solver):
    assert solver is not None


def test_compute_returns_finite_final_velocity(solver, track):
    SS_eval, KK_eval = track
    v_initial = 20.0

    result = solver.compute(SS_eval, KK_eval, v_initial)

    assert math.isfinite(result)
    assert result > 0.0


def test_evaluate_matches_pointwise_eval(solver, track):
    SS_eval, KK_eval = track
    solver.compute(SS_eval, KK_eval, 20.0)

    AX, AY, V = solver.evaluate(SS_eval)

    assert len(AX) == len(AY) == len(V) == len(SS_eval)
    for i, s in enumerate(SS_eval):
        assert V[i] == pytest.approx(solver.evalV(s))
        assert AX[i] == pytest.approx(solver.evalAx(s))
        assert AY[i] == pytest.approx(solver.evalAy(s))


def test_velocity_profile_respects_initial_condition(solver, track):
    SS_eval, KK_eval = track
    v_initial = 20.0
    solver.compute(SS_eval, KK_eval, v_initial)

    assert solver.evalV(SS_eval[0]) == pytest.approx(v_initial, rel=1e-6)


def test_velocity_never_exceeds_vmax(solver, track):
    SS_eval, KK_eval = track
    v_initial = 20.0
    solver.compute(SS_eval, KK_eval, v_initial)

    for s in SS_eval:
        assert solver.evalV(s) <= solver.evalVmax(s) + 1e-6


def test_time_domain_evaluation_is_consistent(solver, track):
    SS_eval, KK_eval = track
    solver.compute_timing(SS_eval, KK_eval, 20.0)

    t0 = 0.0
    s0 = solver.evalS(t0)
    assert s0 == pytest.approx(SS_eval[0], abs=1e-3)
    assert solver.evalV_t(t0) == pytest.approx(solver.evalV(s0), rel=1e-3)
