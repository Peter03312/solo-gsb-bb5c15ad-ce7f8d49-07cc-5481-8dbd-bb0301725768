"""API tests for POST /api/v1/mass-calibrate."""

import math

import pytest
from fastapi.testclient import TestClient

from app.config import settings
from app.calibration import theoretical_mz


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setattr(settings, "db_path", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


ION_SPECS = [
    ("CH4O", "+H", 1),
    ("C2H6O", "+H", 1),
    ("C6H12O6", "+H", 1),
    ("C2H4O2", "-H", -1),
    ("C5H5N5", "+H", 1),
]

SLOPE = 1.000003
INTERCEPT = -0.0006


def _point(spec, slope=SLOPE, intercept=INTERCEPT, factor=1.0):
    formula, adduct, charge = spec
    observed = (slope * theoretical_mz(formula, adduct, charge) + intercept) * factor
    return {
        "observed_mz": observed,
        "formula": formula,
        "adduct": adduct,
        "charge": charge,
    }


def _payload(specs=ION_SPECS, mz_values=None, tolerance_ppm=5.0, **point_overrides):
    points = []
    for i, spec in enumerate(specs):
        point = _point(spec, **point_overrides)
        points.append(point)
    if mz_values is None:
        mz_values = [theoretical_mz(*specs[0]) * SLOPE + INTERCEPT]
    return {
        "calibration_points": points,
        "tolerance_ppm": tolerance_ppm,
        "mz": mz_values,
    }


def test_basic_calibration_recovers_model(client):
    payload = _payload(
        mz_values=[181.07, 91.039, 179.056],
    )
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["slope"] == pytest.approx(SLOPE, abs=1e-8)
    assert body["intercept"] == pytest.approx(INTERCEPT, abs=1e-7)
    assert body["calibration_point_count"] == 5
    assert body["valid_point_count"] == 5
    assert body["inlier_count"] == 5
    assert len(body["calibration_points"]) == 5
    # Request order preserved.
    assert [p["formula"] for p in body["calibration_points"]] == [
        s[0] for s in ION_SPECS
    ]
    for point in body["calibration_points"]:
        assert point["valid"] is True
        assert point["inlier"] is True
        assert point["theoretical_mz"] > 0
        assert abs(point["residual_ppm"]) < 0.1
        assert point["error"] is None


def test_calibrated_mz_inverse_transform_and_order(client):
    true_values = [50.0, 181.070665, 500.0]
    payload = _payload(
        mz_values=[s * SLOPE + INTERCEPT for s in true_values],
    )
    body = client.post("/api/v1/mass-calibrate", json=payload).json()
    assert len(body["calibrated_mz"]) == 3
    for item, true in zip(body["calibrated_mz"], true_values):
        assert item["calibrated_mz"] == pytest.approx(true, abs=1e-4)
    # Original m/z echoed and ordering retained.
    submitted = [s * SLOPE + INTERCEPT for s in true_values]
    assert [item["mz"] for item in body["calibrated_mz"]] == pytest.approx(submitted)


def test_mzs_alias_accepted(client):
    payload = _payload()
    payload["mzs"] = payload.pop("mz")
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 200, resp.text
    assert "calibrated_mz" in resp.json()


def test_mz_values_alias_accepted(client):
    payload = _payload()
    payload["mz_values"] = payload.pop("mz")
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 200, resp.text


def test_invalid_point_is_part_of_200_response(client):
    payload = _payload(specs=ION_SPECS[:4])
    payload["calibration_points"].append(
        {
            "observed_mz": 500.0,
            "formula": "C6Xx",
            "adduct": "+H",
            "charge": 1,
        }
    )
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["calibration_point_count"] == 5
    assert body["valid_point_count"] == 4
    bad = body["calibration_points"][4]
    assert bad["index"] == 4
    assert bad["valid"] is False
    assert bad["inlier"] is None
    assert bad["theoretical_mz"] is None
    assert bad["error"] == "invalid_formula"
    assert "Xx" in bad["message"]


def test_outlier_marked_and_excluded(client):
    specs = ION_SPECS[:5]
    points = [_point(spec) for spec in specs]
    # Add a strong outlier as the last point.
    points.append(
        {
            "observed_mz": theoretical_mz("H2O", "+H", 1) * 0.9,
            "formula": "H2O",
            "adduct": "+H",
            "charge": 1,
        }
    )
    payload = {
        "calibration_points": points,
        "tolerance_ppm": 10.0,
        "mz": [50.0 * SLOPE + INTERCEPT],
    }
    body = client.post("/api/v1/mass-calibrate", json=payload).json()
    assert body["inlier_count"] == 5
    assert body["calibration_points"][-1]["inlier"] is False
    assert body["slope"] == pytest.approx(SLOPE, abs=1e-8)


def test_all_points_invalid_returns_422(client):
    payload = {
        "calibration_points": [
            {"observed_mz": 1.0, "formula": "C6Xx", "adduct": "+H", "charge": 1},
            {"observed_mz": 1.0, "formula": "Zz", "adduct": "+H", "charge": 1},
            {"observed_mz": 1.0, "formula": "H2O", "adduct": "-99H", "charge": 1},
        ],
        "tolerance_ppm": 5.0,
        "mz": [100.0],
    }
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "no_valid_points"


def test_too_few_inliers_returns_422(client):
    specs = ION_SPECS[:3]
    points = []
    for i, spec in enumerate(specs):
        x = theoretical_mz(*spec)
        mz = x if i == 0 else x * 1.01
        points.append(
            {
                "observed_mz": mz,
                "formula": spec[0],
                "adduct": spec[1],
                "charge": spec[2],
            }
        )
    payload = {
        "calibration_points": points,
        "tolerance_ppm": 1.0,
        "mz": [100.0],
    }
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "too_few_inliers"


def test_non_positive_slope_returns_422(client):
    specs = ION_SPECS[:3]
    points = [
        {
            "observed_mz": -2.0 * theoretical_mz(*spec) + 600.0,
            "formula": spec[0],
            "adduct": spec[1],
            "charge": spec[2],
        }
        for spec in specs
    ]
    payload = {
        "calibration_points": points,
        "tolerance_ppm": 1e7,
        "mz": [100.0],
    }
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "non_positive_slope"


def test_point_count_bounds(client):
    base = {
        "tolerance_ppm": 5.0,
        "mz": [100.0],
    }
    too_few = dict(
        base,
        calibration_points=[
            {
                "observed_mz": theoretical_mz(*spec),
                "formula": spec[0],
                "adduct": spec[1],
                "charge": spec[2],
            }
            for spec in ION_SPECS[:2]
        ],
    )
    assert client.post("/api/v1/mass-calibrate", json=too_few).status_code == 422

    spec = ION_SPECS[0]
    too_many = dict(
        base,
        calibration_points=[
            {
                "observed_mz": theoretical_mz(*spec),
                "formula": spec[0],
                "adduct": spec[1],
                "charge": spec[2],
            }
            for _ in range(201)
        ],
    )
    assert client.post("/api/v1/mass-calibrate", json=too_many).status_code == 422


def test_charge_validation(client):
    point = _point(ION_SPECS[0])
    point["charge"] = 0
    payload = {
        "calibration_points": [
            point,
            _point(ION_SPECS[1]),
            _point(ION_SPECS[2]),
        ],
        "tolerance_ppm": 5.0,
        "mz": [100.0],
    }
    assert client.post("/api/v1/mass-calibrate", json=payload).status_code == 422


def test_tolerance_and_mz_validation(client):
    payload = _payload()
    payload["tolerance_ppm"] = 0.0
    assert client.post("/api/v1/mass-calibrate", json=payload).status_code == 422

    payload = _payload()
    payload["mz"] = []
    assert client.post("/api/v1/mass-calibrate", json=payload).status_code == 422

    payload = _payload()
    payload["mz"] = [-1.0]
    assert client.post("/api/v1/mass-calibrate", json=payload).status_code == 422

    payload = _payload()
    payload.pop("mz")
    assert client.post("/api/v1/mass-calibrate", json=payload).status_code == 422


def test_observed_mz_must_be_positive(client):
    point = _point(ION_SPECS[0])
    point["observed_mz"] = 0.0
    payload = {
        "calibration_points": [
            point,
            _point(ION_SPECS[1]),
            _point(ION_SPECS[2]),
        ],
        "tolerance_ppm": 5.0,
        "mz": [100.0],
    }
    assert client.post("/api/v1/mass-calibrate", json=payload).status_code == 422


def test_large_legal_observations_fit_and_calibrate(client):
    # Regression: astronomical but finite observed m/z (~1e154) used to
    # overflow the inlier least-squares accumulators (inf slope, NaN
    # residuals and a spurious 422). The representable model must come back.
    base = 10**153
    points = []
    for mult in (1, 2, 3):
        formula = f"C{base * mult}"
        x = theoretical_mz(formula, "+H", 1)
        points.append(
            {
                "observed_mz": x,
                "formula": formula,
                "adduct": "+H",
                "charge": 1,
            }
        )
    payload = {
        "calibration_points": points,
        "tolerance_ppm": 5.0,
        "mz": [theoretical_mz(f"C{base}", "+H", 1)],
    }
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert math.isfinite(body["slope"]) and body["slope"] > 0
    assert math.isfinite(body["intercept"])
    assert body["inlier_count"] == 3
    for point in body["calibration_points"]:
        assert point["valid"] is True
        assert point["inlier"] is True
        assert math.isfinite(point["theoretical_mz"])
        assert math.isfinite(point["residual_mz"])
        assert math.isfinite(point["residual_ppm"])
    item = body["calibrated_mz"][0]
    assert math.isfinite(item["calibrated_mz"])
    assert item["calibrated_mz"] == pytest.approx(
        theoretical_mz(f"C{base}", "+H", 1), rel=1e-6
    )


def test_large_negative_intercept_calibration_is_returned(client):
    # Regression: a finite fitted model with a large negative intercept made
    # the inverse subtraction overflow before the division, so an otherwise
    # representable calibrated m/z came back as Infinity.
    slope, intercept = 2.0, -1.7e308
    xtarget = 0.95 * 1.7976931348623157e308 / slope
    c0 = int(xtarget / 12.0)
    points = []
    xs = []
    for k in range(3):
        formula = f"C{c0 + k * 10**292}"
        x = theoretical_mz(formula, "+H", 1)
        y = slope * x + intercept
        assert math.isfinite(y) and y > 0
        xs.append(x)
        points.append(
            {"observed_mz": y, "formula": formula, "adduct": "+H", "charge": 1}
        )
    payload = {
        "calibration_points": points,
        "tolerance_ppm": 1e4,
        "mz": [1.7e308],
    }
    resp = client.post("/api/v1/mass-calibrate", json=payload)
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert math.isfinite(body["slope"])
    assert math.isfinite(body["intercept"])
    calibrated = body["calibrated_mz"][0]["calibrated_mz"]
    assert math.isfinite(calibrated)
    assert calibrated == pytest.approx(1.7e308, rel=1e-6)


def test_request_logged(client):
    from app import db

    payload = _payload()
    assert client.post("/api/v1/mass-calibrate", json=payload).status_code == 200
    with db._require_conn() as conn:
        rows = conn.execute(
            "SELECT formula, ok, message FROM request_log WHERE formula LIKE "
            "'mass-calibrate%'"
        ).fetchall()
    assert len(rows) == 1
    assert rows[0][1] == 1
    assert "inliers" in rows[0][2]
