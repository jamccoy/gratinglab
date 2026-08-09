"""GUI form parsing.

Runs headlessly -- no display, no tkinter -- because this module deliberately
holds every piece of GUI logic that could be wrong in an interesting way. The
widget layer gets none.
"""

import numpy as np
import pytest

from gratinglab.gui.scalar_options import ScalarOptionsState, build_options
from gratinglab.gui.state import (
    ANGLE_LABELS,
    MOUNTS,
    PROFILE_KINDS,
    FormErrors,
    FormState,
    build,
    validate,
)
from gratinglab.illumination import Illumination
from gratinglab.io.ggp import write_ggp
from gratinglab.profiles import Blazed, FromProfileData, Lamellar, Sinusoidal
from gratinglab.solvers import scalar


class TestDefaults:
    def test_the_default_form_is_valid(self):
        """A user who opens the app and presses Solve must get a result."""
        assert validate(FormState()) == ()

    def test_the_default_is_the_soft_xray_off_plane_case(self):
        parsed = build(FormState())
        assert parsed.problem.period == pytest.approx(315.15)
        assert parsed.illumination.gamma_deg == pytest.approx(1.5)
        assert not parsed.illumination.is_in_plane

    def test_the_default_actually_solves(self):
        """End to end through the real solver, not a mock.

        `build()` no longer carries solver options -- see
        `gratinglab.gui.scalar_options` -- so scalar's own defaults are built
        separately, the way a solver tab does.
        """
        parsed = build(FormState())
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(),
        )
        scan = scalar.solve(
            parsed.problem,
            parsed.illumination,
            parsed.wavelengths,
            **options,
        )
        assert len(scan) == 200
        assert np.isfinite(scan.efficiency).all()


class TestProfiles:
    def test_blazed(self):
        profile = build(FormState(profile_kind="Blazed")).problem.profile
        assert isinstance(profile, Blazed)
        assert profile.blaze_angle == pytest.approx(29.5)

    def test_lamellar(self):
        parsed = build(
            FormState(profile_kind="Lamellar", depth_fraction="0.2", duty_cycle="0.4")
        )
        profile = parsed.problem.profile
        assert isinstance(profile, Lamellar)
        assert profile.duty_cycle == pytest.approx(0.4)

    def test_sinusoidal(self):
        profile = build(
            FormState(profile_kind="Sinusoidal", depth_fraction="0.15")
        ).problem.profile
        assert isinstance(profile, Sinusoidal)

    def test_from_file(self, tmp_path):
        path = write_ggp(tmp_path / "p.ggp", t=[0.0, 0.4, 1.0], y=[0.0, 0.25, 0.0])
        profile = build(
            FormState(profile_kind="From file", profile_path=str(path))
        ).problem.profile
        assert isinstance(profile, FromProfileData)
        assert profile.depth == pytest.approx(0.25, abs=1e-5)

    def test_missing_file_is_a_field_error_not_a_traceback(self):
        errors = validate(
            FormState(profile_kind="From file", profile_path="/nope/missing.ggp")
        )
        assert [e.field for e in errors] == ["profile_path"]
        assert "no such file" in errors[0].message

    def test_empty_path_asks_for_one(self):
        errors = validate(FormState(profile_kind="From file", profile_path=""))
        assert "choose a profile file" in errors[0].message

    def test_a_corrupt_file_reports_the_readers_message(self, tmp_path):
        bad = tmp_path / "bad.ggp"
        bad.write_text("3 0 - Polygonal type\nPeriod: 1 PSC: 1\n0.0 0.0\n0.5\n1.0 0.0\n")
        errors = validate(FormState(profile_kind="From file", profile_path=str(bad)))
        assert errors and errors[0].field == "profile_path"
        assert "x y" in errors[0].message

    def test_unknown_kind(self):
        errors = validate(FormState(profile_kind="Helical"))
        assert any(e.field == "profile_kind" for e in errors)


class TestMounts:
    def test_classical_is_in_plane(self):
        ill = build(FormState(mount="Classical", alpha="10")).illumination
        assert ill.is_in_plane
        assert ill.alpha_deg == pytest.approx(10.0)

    def test_offplane_maps_gamma_to_graze(self):
        ill = build(
            FormState(mount="Off-plane", alpha="25", gamma="1.5")
        ).illumination
        assert ill.gamma_deg == pytest.approx(1.5)
        assert ill.alpha_deg == pytest.approx(25.0)

    def test_conical_maps_the_two_angles_to_theta_and_phi(self):
        ill = build(FormState(mount="Conical", alpha="30", gamma="0")).illumination
        assert ill == Illumination.conical(
            theta=30.0, phi=0.0, polarization="unpolarized"
        )

    def test_every_mount_has_angle_labels(self):
        """The UI relabels the same two widgets; a missing entry would crash it."""
        assert set(ANGLE_LABELS) == set(MOUNTS)

    def test_classical_hides_the_second_angle(self):
        assert ANGLE_LABELS["Classical"][1] is None

    def test_unknown_mount(self):
        assert any(e.field == "mount" for e in validate(FormState(mount="Littrow")))


class TestNumericValidation:
    @pytest.mark.parametrize(
        "field,value",
        [
            ("period", "0"),
            ("period", "-5"),
            ("period", "abc"),
            ("period", ""),
            ("blaze_angle", "90"),
            ("blaze_angle", "0"),
            ("wavelength_start", "0"),
            ("wavelength_count", "0"),
        ],
    )
    def test_rejects_out_of_range_or_unparseable(self, field, value):
        errors = validate(FormState().with_field(field, value))
        assert any(e.field == field for e in errors), f"{field}={value} accepted"

    def test_rejects_non_integer_counts(self):
        errors = validate(FormState(wavelength_count="10.5"))
        assert any("whole number" in e.message for e in errors)

    def test_rejects_nan_and_inf(self):
        for value in ("nan", "inf", "-inf"):
            errors = validate(FormState(period=value))
            assert errors, f"{value} accepted as a period"

    def test_rejects_a_backwards_wavelength_range(self):
        errors = validate(FormState(wavelength_start="5", wavelength_stop="1"))
        assert any("must exceed the start" in e.message for e in errors)

    def test_accepts_surrounding_whitespace(self):
        assert build(FormState(period="  315.15  ")).problem.period == pytest.approx(
            315.15
        )

    def test_reports_every_bad_field_at_once(self):
        """Fixing a form one error per attempt is a miserable experience."""
        errors = validate(
            FormState(period="oops", blaze_angle="-1", wavelength_count="0")
        )
        assert {e.field for e in errors} >= {
            "period",
            "blaze_angle",
            "wavelength_count",
        }

    def test_error_messages_name_the_bound(self):
        message = validate(FormState(blaze_angle="95"))[0].message
        assert "90" in message


class TestFormState:
    def test_with_field_does_not_mutate(self):
        original = FormState()
        changed = original.with_field("period", "160")
        assert original.period == "315.15"
        assert changed.period == "160"

    def test_every_offered_profile_kind_builds(self, tmp_path):
        path = write_ggp(tmp_path / "p.ggp", t=[0.0, 0.4, 1.0], y=[0.0, 0.25, 0.0])
        for kind in PROFILE_KINDS:
            form = FormState(profile_kind=kind, profile_path=str(path))
            assert validate(form) == (), f"{kind} does not build with defaults"

    def test_every_offered_mount_builds(self):
        for mount in MOUNTS:
            form = FormState(mount=mount, alpha="20", gamma="30")
            assert validate(form) == (), f"{mount} does not build"

    def test_lambda_over_period_is_reported(self):
        parsed = build(FormState())
        assert parsed.lambda_over_period == pytest.approx(5.0 / 315.15)


class TestNoWidgetDependency:
    def test_state_module_imports_no_toolkit(self):
        """The whole reason this module exists separately.

        Checked against the parsed imports rather than the source text -- the
        word appears in the module docstring, which is not a dependency.
        """
        import ast
        import inspect

        import gratinglab.gui.state as state

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(state))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert "tkinter" not in imported
        assert not {"matplotlib", "PySide6", "PyQt6"} & imported

    def test_state_is_importable_without_a_display(self):
        """CI has no display; this module must still import and run there."""
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-c",
             "import gratinglab.gui.state as s; "
             "assert s.validate(s.FormState()) == (); print('ok')"],
            capture_output=True, text=True, env={"PATH": "/usr/bin:/bin"},
        )
        assert result.returncode == 0, result.stderr
        assert "ok" in result.stdout

    def test_formerrors_carries_structured_errors(self):
        with pytest.raises(FormErrors) as excinfo:
            build(FormState(period="bad"))
        assert excinfo.value.errors[0].field == "period"
        assert "period" in str(excinfo.value)
