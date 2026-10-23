"""
Profile flattening, the stage that does move the answer.

`level_grooves` fits its background through detected features, which makes the
size of that feature set load-bearing. Two ways it used to go wrong silently:

- Edge exclusion ran *after* detection and overwrote the feature arrays, so a
  trace whose features all sat near the ends was fitted on whatever remained -
  sometimes a single point - while the code printed "using all".
- The guard was `len(feature_x) < 2`, but `poly_order` defaults to 2 and a
  degree-k polynomial needs k+1 points. Two features and a quadratic gave a
  rank-deficient background behind a numpy RankWarning.

Both produce a background that looks like a fit and is not one, and because
flattening sets the plane every blaze angle on the scan is quoted against, a
bad background is a bias on every groove rather than noise on one. These tests
pin the handling.

Each fallback raises a `UserWarning`. The suite runs with
`filterwarnings = error`, so a test that provokes one has to say so, and any
other test that starts flattening by fallback fails rather than passing
quietly.

The profiles here are built rather than measured: a few smooth bumps at known
positions, so the feature set is exactly what the test intends and the
expected background can be written down.
"""
import numpy as np
import pytest

from gratinglab.metrology.core.processing import flatten_profile

#: x in micrometres, 1 nm per sample - `flatten_profile` reads dx in nm as
#: `(x[1] - x[0]) * 1000`, and the feature separation gate is in periods.
DX_UM = 0.001
PERIOD_NM = 100.0


def _profile(centres_um, heights=None, n=1001, width_um=0.01):
    """A baseline of zero with a smooth bump at each centre."""
    x = np.arange(n) * DX_UM
    heights = [1.0] * len(centres_um) if heights is None else heights
    y = np.zeros_like(x)
    for centre, height in zip(centres_um, heights):
        y = y + height * np.exp(-(((x - centre) / width_um) ** 2))
    return x, y


def _background_through(x, feature_x, feature_y, order):
    """The polynomial a correct implementation must have fitted."""
    return np.polyval(np.polyfit(feature_x, feature_y, order), x)


class TestEdgeExclusionDoesNotStarveTheFit:
    """
    Excluding the edges is an improvement to the fit, not a licence to break it.

    When exclusion would leave too few features for the requested order, the
    full set is used and the warning says so. The old code said so and then
    fitted the remnant anyway.
    """

    CENTRES = (0.05, 0.5, 0.95)

    def test_all_features_are_used_when_exclusion_would_leave_too_few(self):
        x, y = _profile(self.CENTRES)

        with pytest.warns(UserWarning, match='using all of them'):
            flat, background = flatten_profile(
                x, y, method='level_grooves', poly_order=2,
                exclude_edges=0.2, feature='peaks', period_nm=PERIOD_NM)

        # exclude_edges=0.2 keeps only 0.2 <= x <= 0.8, i.e. the middle bump.
        # A degree-2 background needs three points, so all three are used.
        expected = _background_through(
            x, np.array(self.CENTRES), np.ones(3), 2)

        assert np.allclose(background, expected, atol=1e-9), (
            "the background was not fitted through every detected feature")
        assert np.allclose(flat, y - expected, atol=1e-9)

    def test_the_warning_says_what_actually_happened(self):
        x, y = _profile(self.CENTRES)

        with pytest.warns(UserWarning) as caught:
            flatten_profile(x, y, method='level_grooves', poly_order=2,
                            exclude_edges=0.2, feature='peaks',
                            period_nm=PERIOD_NM)

        message = str(caught[0].message)
        assert 'edge exclusion leaves 1 feature(s)' in message, message
        assert 'using all of them' in message, message

    def test_exclusion_still_bites_when_enough_features_survive(self):
        """
        The fix must not quietly disable the setting it guards.

        Eleven evenly spaced bumps, the outer ones raised: excluding the edges
        has to change the background, or edge exclusion has stopped working.
        """
        centres = [0.05 + 0.09 * i for i in range(11)]
        heights = [2.0] + [1.0] * 9 + [2.0]
        x, y = _profile(centres, heights)

        # Neither call may warn: eleven features survive either way, so a
        # warning here would mean the fit was starved after all.
        _, excluded = flatten_profile(
            x, y, method='level_grooves', poly_order=2,
            exclude_edges=0.2, feature='peaks', period_nm=PERIOD_NM)
        _, everything = flatten_profile(
            x, y, method='level_grooves', poly_order=2,
            exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)

        assert not np.allclose(excluded, everything, atol=1e-6), (
            "excluding the edges made no difference, so it was not applied")


class TestOrderIsReducedRatherThanRankDeficient:
    """
    Two features cannot determine a quadratic.

    numpy will return something regardless, behind a RankWarning that nothing
    in this pipeline was reading. Reducing the order gives the best background
    the data actually supports.
    """

    def test_two_features_give_the_line_through_them(self):
        x, y = _profile((0.25, 0.75))

        with pytest.warns(UserWarning, match='degree-1 background'):
            _, background = flatten_profile(
                x, y, method='level_grooves', poly_order=2,
                exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)

        expected = _background_through(x, np.array((0.25, 0.75)), np.ones(2), 1)
        assert np.allclose(background, expected, atol=1e-9)

        # A straight line, not a curve that happened to pass through two points.
        assert np.allclose(np.diff(background, 2), 0.0, atol=1e-12)

    def test_the_reduction_is_reported(self):
        x, y = _profile((0.25, 0.75))

        with pytest.warns(UserWarning) as caught:
            flatten_profile(x, y, method='level_grooves', poly_order=2,
                            exclude_edges=0.0, feature='peaks',
                            period_nm=PERIOD_NM)

        message = str(caught[0].message)
        assert 'support only a degree-1 background' in message, message


class TestTheUnaffectedPathsStayPut:
    """
    The fix touches one branch. These pin the behaviour around it.
    """

    def test_ample_features_fit_the_requested_order(self):
        centres = [0.05 + 0.09 * i for i in range(11)]
        x, y = _profile(centres)

        _, background = flatten_profile(
            x, y, method='level_grooves', poly_order=2,
            exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)

        expected = _background_through(
            x, np.array(centres), np.ones(len(centres)), 2)
        assert np.allclose(background, expected, atol=1e-9)

    def test_too_few_features_still_falls_back_to_linear(self):
        """
        One bump cannot be levelled on its features, and the outer guard
        catches that before the branch this change touches.
        """
        x, y = _profile((0.5,))

        with pytest.warns(UserWarning, match='using linear'):
            _, background = flatten_profile(
                x, y, method='level_grooves', poly_order=2,
                exclude_edges=0.0, feature='peaks', period_nm=PERIOD_NM)
        assert np.allclose(background, np.polyval(np.polyfit(x, y, 1), x))

    def test_level_grooves_still_demands_a_period(self):
        x, y = _profile((0.25, 0.75))
        with pytest.raises(ValueError, match='period_nm'):
            flatten_profile(x, y, method='level_grooves', period_nm=None)
