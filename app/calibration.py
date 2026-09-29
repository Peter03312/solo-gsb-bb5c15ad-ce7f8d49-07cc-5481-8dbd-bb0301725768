"""Mass-axis calibration from reference calibration points.

Each calibration point pairs an observed monoisotopic m/z with the neutral
formula, adduct string and signed charge of the reference ion; its
theoretical monoisotopic m/z is computed with the same formula parsing,
adduct application and electron-mass rule as ``/api/v1/ion-envelope``::

    m/z = (monoisotopic_mass(adjusted) - charge * m_e) / |charge|

The calibration model is an affine map from theoretical m/z (x) to observed
m/z (y): ``y = slope * x + intercept``. A robust initial model is built with
the Theil–Sen/pairwise construction: the initial slope is the median of the
slopes defined by every pair of valid points with distinct theoretical m/z,
and the initial intercept is the median of that pair's intercept (the same
line passes through both points, so the intercept computed from either one
is identical). Points whose initial-model residual is within
``tolerance_ppm`` (``|residual / theoretical| * 1e6 <= tolerance_ppm``) are
the inliers; a single ordinary least-squares line is then fitted to the
inliers and reported as the final model. Calibrating an observed m/z applies
the inverse transform ``(mz - intercept) / slope``.

The fit aborts with :class:`CalibrationError` when no point is valid, when
no pair with distinct theoretical m/z exists, when fewer than three points
are inliers, or when the least-squares slope is not positive (the inverse
transform needs a strictly positive slope).

Observed m/z values may be as large as the float range allows. Every
floating-point step of the fit (pairwise slope/intercept, median averages,
least-squares accumulation, residual evaluation and the inverse transform)
is therefore evaluated with exponent-scaled arithmetic so that an
intermediate term which merely passes through the overflow regime cannot
produce a spurious infinity or NaN: only results that are genuinely outside
the finite float range overflow.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .adduct import AdductError, apply_adduct, parse_adduct
from .envelope import ELECTRON_MASS_DA, monoisotopic_mass
from .parser import FormulaError, parse_formula

# A residual exactly on the tolerance boundary is an inlier; this guard only
# absorbs ordinary floating-point round-off at that boundary.
_PPM_EPS = 1e-9

MIN_POINTS = 3
MAX_POINTS = 200
MIN_INLIERS = 3


class CalibrationError(ValueError):
    """Calibration could not be built from the submitted points."""

    def __init__(self, error: str, message: str) -> None:
        super().__init__(message)
        self.error = error
        self.message = message


@dataclass(frozen=True)
class CalibrationPoint:
    """One reference calibration point in request order."""

    index: int
    observed_mz: float
    formula: str
    adduct: str
    charge: int


@dataclass
class PointResidual:
    """Per-point outcome: theoretical m/z, final-model residual and the
    inlier/outlier flag assigned by the robust initial model."""

    index: int
    valid: bool
    observed_mz: float
    formula: str
    adduct: str
    charge: int
    theoretical_mz: float | None = None
    residual_mz: float | None = None
    residual_ppm: float | None = None
    inlier: bool | None = None
    error: str | None = None
    message: str | None = None


@dataclass
class MassCalibrationFit:
    """Final affine calibration model plus per-point diagnostics."""

    slope: float
    intercept: float
    valid_point_count: int
    inlier_count: int
    points: list[PointResidual]


def theoretical_mz(formula: str, adduct: str, charge: int) -> float:
    """Theoretical monoisotopic m/z of a reference ion.

    Reuses the existing formula parser, adduct application and the
    ``(mass - charge * m_e) / |charge|`` electron-mass rule shared with the
    ion-envelope endpoints. Raises :class:`FormulaError` or
    :class:`AdductError` exactly as those parsers do.
    """
    parsed = parse_formula(formula)
    pieces = parse_adduct(adduct)
    adjusted = apply_adduct(parsed, pieces)
    return (monoisotopic_mass(adjusted) - charge * ELECTRON_MASS_DA) / abs(charge)


def _stable_fma(a: float, b: float, c: float) -> float:
    """Fused-multiply-add ``a*b + c`` without spurious overflow.

    A plain ``a * b + c`` can overflow even when the true sum is finite and
    representable (e.g. a huge product that is nearly cancelled by ``c``).
    The operands are split into a common exponent scale so every partial
    term stays in the finite range; the scale is restored afterwards. An
    infinity is returned only when the mathematical result itself exceeds
    the finite float range. Inputs that are already non-finite fall through
    to the ordinary IEEE expression.
    """
    if not (math.isfinite(a) and math.isfinite(b) and math.isfinite(c)):
        return a * b + c
    if a == 0.0 or b == 0.0:
        return c + 0.0  # signed zero of the addend, never an overflow

    pa, ea = math.frexp(a)
    pb, eb = math.frexp(b)
    pc, ec = math.frexp(c) if c != 0.0 else (0.0, 0)

    # Work at the scale of the product: a*b = pa*pb * 2^(ea+eb) and
    # c = pc * 2^ec, and both are evaluated as fractions of 2^q.
    q = max(ea + eb, ec)
    scaled_product = math.ldexp(pa * pb, ea + eb - q)
    scaled_addend = (
        math.ldexp(pc, ec - q) if c != 0.0 else 0.0
    )
    try:
        return math.ldexp(scaled_product + scaled_addend, q)
    except OverflowError:
        return math.copysign(math.inf, a * b + c)


def _stable_sum2(a: float, b: float) -> float:
    """Stable two-operand difference used for residual evaluation.

    Residuals are differences of similarly sized quantities. For finite
    operands the IEEE ``+`` (with the second operand negated) already
    overflows only when the mathematical result genuinely exceeds the
    finite range: operands of opposite sign can never overflow, since the
    magnitude of their sum is no larger than the larger operand. Keeping
    this as an explicit call site documents the invariant and lets
    non-finite operands follow ordinary IEEE propagation.
    """
    return a + b


def _median(values: list[float]) -> float:
    """Median of a non-empty list (average of the two middle values on ties).

    The even-length average is taken as the mean of two half-size terms so
    that averaging two huge, near-cancelling middle values cannot overflow.
    """
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return _stable_sum2(ordered[mid - 1] / 2.0, ordered[mid] / 2.0)


def _point_result_invalid(point: CalibrationPoint, exc: Exception) -> PointResidual:
    if isinstance(exc, FormulaError):
        error = "invalid_formula"
    else:
        error = "invalid_adduct"
    return PointResidual(
        index=point.index,
        valid=False,
        observed_mz=point.observed_mz,
        formula=point.formula,
        adduct=point.adduct,
        charge=point.charge,
        error=error,
        message=str(exc.args[0]) if exc.args else str(exc),
    )


def calibrate_mass_axis(
    calibration_points: list[CalibrationPoint],
    tolerance_ppm: float,
) -> MassCalibrationFit:
    """Build the robust-then-least-squares mass-axis calibration model.

    Invalid points (formula/adduct that fail to parse or apply) are reported
    in place and excluded from the fit. The fit raises
    :class:`CalibrationError` for the failure cases described in the module
    docstring.
    """
    points: list[PointResidual] = []
    valid: list[tuple[float, float]] = []  # (theoretical x, observed y)

    for point in calibration_points:
        try:
            x = theoretical_mz(point.formula, point.adduct, point.charge)
        except (FormulaError, AdductError) as exc:
            points.append(_point_result_invalid(point, exc))
            continue
        points.append(
            PointResidual(
                index=point.index,
                valid=True,
                observed_mz=point.observed_mz,
                formula=point.formula,
                adduct=point.adduct,
                charge=point.charge,
                theoretical_mz=x,
            )
        )
        valid.append((x, point.observed_mz))

    if not valid:
        raise CalibrationError(
            "no_valid_points",
            "No valid calibration points: every point failed formula or "
            "adduct parsing; see the per-point error/message fields",
        )

    # Robust initial model: median of pairwise slopes (Theil–Sen) and the
    # median intercept of the line through each pair. Only pairs with
    # distinct theoretical m/z define a finite slope. A pair whose slope or
    # intercept is non-finite (e.g. a theoretical m/z beyond the float
    # range) cannot define a usable affine line and is skipped.
    slopes: list[float] = []
    intercepts: list[float] = []
    for i in range(len(valid) - 1):
        x_i, y_i = valid[i]
        for j in range(i + 1, len(valid)):
            x_j, y_j = valid[j]
            if x_i == x_j or not (math.isfinite(x_i) and math.isfinite(x_j)):
                continue
            pair_slope = (y_j - y_i) / (x_j - x_i)
            # y_i - pair_slope * x_i evaluated stably: for very close
            # (x_i, x_j) the product can run through the overflow regime
            # even though the true intercept stays in range.
            pair_intercept = _stable_fma(-pair_slope, x_i, y_i)
            if not (math.isfinite(pair_slope) and math.isfinite(pair_intercept)):
                continue
            slopes.append(pair_slope)
            intercepts.append(pair_intercept)

    if not slopes:
        raise CalibrationError(
            "no_valid_point_pairs",
            "Cannot build an initial calibration model: all valid points "
            "share the same theoretical m/z",
        )

    initial_slope = _median(slopes)
    initial_intercept = _median(intercepts)

    inlier_flags: list[bool] = []
    for x, y in valid:
        # Stabilized prediction: initial_slope * x + initial_intercept can
        # otherwise overflow for large theoretical m/z even when the
        # subtraction against y leaves a small residual.
        predicted = _stable_fma(initial_slope, x, initial_intercept)
        residual = _stable_sum2(y, -predicted)
        if math.isfinite(residual) and x > 0.0:
            residual_ppm = residual / x * 1e6
        else:
            residual_ppm = math.inf
        inlier_flags.append(abs(residual_ppm) <= tolerance_ppm + _PPM_EPS)

    inlier_count = sum(inlier_flags)
    if inlier_count < MIN_INLIERS:
        raise CalibrationError(
            "too_few_inliers",
            f"Only {inlier_count} of {len(valid)} valid calibration points "
            f"are inliers within {tolerance_ppm:g} ppm; at least "
            f"{MIN_INLIERS} are required (widen tolerance_ppm or replace "
            f"outlying points)",
        )

    inlier_x = [x for (x, _), flag in zip(valid, inlier_flags) if flag]
    inlier_y = [y for (_, y), flag in zip(valid, inlier_flags) if flag]

    # Ordinary least squares in a common exponent-scaled space. Naive
    # sum-of-squares accumulation overflows for legitimate large-magnitude
    # points (each (x-mean)**2 term can exceed the float range, and inf-inf
    # cancellation in sxy/sxx would otherwise yield a NaN slope and NaN
    # residuals). Normalising x and y by their own binary exponents keeps
    # every accumulated term below 2, and the scale is restored analytically
    # in the slope and intercept.
    sx = math.ldexp(1.0, math.frexp(max(abs(v) for v in inlier_x))[1] - 1)
    sy = math.ldexp(1.0, math.frexp(max(abs(v) for v in inlier_y))[1] - 1)
    ux = [x / sx for x in inlier_x]
    uy = [y / sy for y in inlier_y]
    if not all(math.isfinite(v) for v in (*ux, *uy)):
        raise CalibrationError(
            "degenerate_inliers",
            "The inlier points contain non-finite theoretical or observed "
            "m/z values, so a least-squares line cannot be fitted",
        )
    n = inlier_count
    mean_u = math.fsum(ux) / n
    mean_v = math.fsum(uy) / n
    suu = math.fsum((u - mean_u) ** 2 for u in ux)
    suv = math.fsum((u - mean_u) * (v - mean_v) for u, v in zip(ux, uy))
    if suu <= 0.0:
        raise CalibrationError(
            "degenerate_inliers",
            "The inlier points share one theoretical m/z, so a "
            "least-squares slope is undefined",
        )
    # slope = (sy/sx) * suv/suu, with the factors combined through ldexp to
    # avoid an intermediate overflow when the true slope is finite.
    factor, factor_exp = math.frexp(sy / sx)
    q, q_exp = math.frexp(suv / suu)
    try:
        slope = math.ldexp(factor * q, factor_exp + q_exp)
    except OverflowError:
        slope = math.copysign(math.inf, suv)
    mean_x = mean_u * sx
    mean_y = mean_v * sy
    intercept = _stable_fma(-slope, mean_x, mean_y)
    if not (math.isfinite(slope) and slope > 0.0):
        raise CalibrationError(
            "non_positive_slope",
            f"The least-squares calibration slope is {slope:.6g}, which is "
            "not a finite positive number; the inverse transform is "
            "undefined (check the calibration points for swapped or "
            "inconsistent m/z values)",
        )
    if not math.isfinite(intercept):
        raise CalibrationError(
            "degenerate_inliers",
            "The least-squares calibration intercept is not finite, so the "
            "fitted model cannot represent a stable calibration",
        )

    # Final-model residuals for every valid point; the inlier/outlier flag
    # stays the one assigned by the robust initial model.
    valid_iter = iter(zip(valid, inlier_flags))
    for result in points:
        if not result.valid:
            continue
        (x, y), flag = next(valid_iter)
        residual = _stable_sum2(y, -_stable_fma(slope, x, intercept))
        result.residual_mz = residual
        if math.isfinite(residual) and math.isfinite(x) and x > 0.0:
            result.residual_ppm = residual / x * 1e6
        else:
            result.residual_ppm = math.inf
        result.inlier = flag

    return MassCalibrationFit(
        slope=slope,
        intercept=intercept,
        valid_point_count=len(valid),
        inlier_count=inlier_count,
        points=points,
    )


def inverse_predict(fit: MassCalibrationFit, mz: float) -> float:
    """Apply the fitted inverse transform to one observed m/z.

    ``(mz - intercept) / slope`` is stabilized so the subtraction cannot
    overflow merely because a negative intercept is being removed: the
    quotient is evaluated in a common binary scale of the three operands.
    The result overflows only when the calibrated m/z itself lies outside
    the finite float range.
    """
    slope = fit.slope
    intercept = fit.intercept
    # Fast path: a subtraction of finite numbers is infinite only if the
    # true difference overflows, in which case a direct quotient is needed.
    diff = _stable_sum2(mz, -intercept)
    if math.isfinite(diff):
        return diff / slope

    if not (math.isfinite(mz) and math.isfinite(intercept) and math.isfinite(slope)):
        return (mz - intercept) / slope

    # Shared scale at the largest operand: all three scaled operands stay
    # below 2, so the subtraction stays finite. The scale cancels in the
    # quotient ((mz-b)/scale)/(slope/scale) == (mz-b)/slope, so the result
    # overflows only when the calibrated m/z itself is out of finite range.
    emax = max(
        math.frexp(mz)[1],
        math.frexp(intercept)[1],
        math.frexp(slope)[1],
    )
    num = math.ldexp(mz, 1 - emax) - math.ldexp(intercept, 1 - emax)
    den = math.ldexp(slope, 1 - emax)
    if den == 0.0:
        if num == 0.0:
            return 0.0
        positive = (num > 0.0) == (slope > 0.0)
        return math.inf if positive else -math.inf
    return num / den


__all__ = [
    "MAX_POINTS",
    "MIN_INLIERS",
    "MIN_POINTS",
    "CalibrationError",
    "CalibrationPoint",
    "MassCalibrationFit",
    "PointResidual",
    "calibrate_mass_axis",
    "inverse_predict",
    "theoretical_mz",
]
