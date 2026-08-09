"""The grating geometry: period, profile, mount, wavelengths.

Extracted from `main_window.py`'s original `_build_inputs`, minus scalar's own
"Scalar solver" groupbox and its Solve/Cancel/Export controls -- those are
solver-specific and stay wherever a solver's own tab lives (`ScalarTab`,
eventually alongside whatever RCWA's tab turns out to need). This panel is
shared: it is built once, not once per solver, because a `Problem` and an
`Illumination` mean the same thing regardless of which solver is about to
consume them.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..state import ANGLE_LABELS, MOUNTS, PROFILE_FIELDS, PROFILE_KINDS, FormState

__all__ = ["GeometryPanel"]


class GeometryPanel(QWidget):
    """Grating, mount and wavelength inputs -- the persistent left sidebar.

    Not wrapped in its own scroll area: how this panel is embedded (today,
    stacked above a solver's own controls under one shared scroller) is the
    caller's decision, not this widget's.
    """

    #: Pressing Enter in a field, or finishing a "Load profile…" pick. This
    #: panel does not know what a solve is, only that something changed
    #: enough that one might be wanted; `MainWindow` connects it to `solve()`.
    solve_requested = Signal()

    #: The form was edited *by a user*. Deliberately distinct from
    #: `solve_requested`: Enter still means solve, and this means only that
    #: the geometry drawing is now stale. `MainWindow` debounces it and
    #: redraws; nothing here reaches a solver.
    #:
    #: Wired to `textEdited`/`activated`, never `textChanged`/
    #: `currentTextChanged`, because those also fire on programmatic
    #: `setText`/`setCurrentText` -- which is how every test and
    #: `_on_mount_change` set fields. Only a real edit should schedule work.
    changed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._fields: dict[str, QWidget] = {}
        self._build()
        self._on_mount_change()
        self._on_profile_change()
        self._check_fields_match_formstate()

    def _build(self) -> None:
        defaults = FormState()
        column = QVBoxLayout(self)
        column.setContentsMargins(0, 0, 0, 0)

        def entry(form: QFormLayout, label: str, key: str) -> QLineEdit:
            widget = QLineEdit(getattr(defaults, key))
            widget.returnPressed.connect(self.solve_requested)
            widget.textEdited.connect(self.changed)
            self._fields[key] = widget
            form.addRow(label, widget)
            return widget

        def combo(form: QFormLayout, label: str, key: str, values, on_change) -> QComboBox:
            widget = QComboBox()
            widget.addItems(list(values))
            widget.setCurrentText(getattr(defaults, key))
            widget.currentTextChanged.connect(lambda _t: on_change())
            widget.activated.connect(lambda _i: self.changed.emit())
            self._fields[key] = widget
            form.addRow(label, widget)
            return widget

        grating = QGroupBox("Grating")
        grating_form = self._grating_form = QFormLayout(grating)
        entry(grating_form, "Period (nm)", "period")
        combo(grating_form, "Profile", "profile_kind", PROFILE_KINDS,
              self._on_profile_change)

        self._profile_rows: dict[str, QWidget] = {}
        for key, label in (
            ("blaze_angle", "Blaze δ (deg)"),
            ("antiblaze_angle", "Anti-blaze (deg)"),
            ("depth_fraction", "Depth / period"),
            ("duty_cycle", "Duty cycle"),
        ):
            self._profile_rows[key] = entry(grating_form, label, key)

        # Not a visible row: the path is set by the file dialog, and shown by
        # the label beneath it. It still has to be a field, because
        # read_form() builds FormState from exactly this dict.
        path_field = QLineEdit("")
        path_field.hide()
        self._fields["profile_path"] = path_field

        self._path_label = QLabel("")
        self._path_label.setWordWrap(True)
        self._path_label.setStyleSheet("color: #555")
        self._load_button = QPushButton("Load profile…")
        self._load_button.clicked.connect(self.load_profile)
        grating_form.addRow(self._load_button)
        grating_form.addRow(self._path_label)

        mount = QGroupBox("Mount")
        mount_form = QFormLayout(mount)
        combo(mount_form, "Geometry", "mount", MOUNTS, self._on_mount_change)
        self._angle_labels: dict[str, QLabel] = {}
        for key in ("alpha", "gamma"):
            widget = QLineEdit(getattr(defaults, key))
            widget.returnPressed.connect(self.solve_requested)
            widget.textEdited.connect(self.changed)
            self._fields[key] = widget
            label = QLabel("")
            mount_form.addRow(label, widget)
            self._angle_labels[key] = label

        surface = QGroupBox("Surface")
        surface_form = QFormLayout(surface)
        # A dropdown, not a text field: the set of materials is finite and
        # known, and a free-form string is what let `coating` mean nothing for
        # four milestones. Empty entry first, because no coating is the
        # ordinary default and gives relative efficiency -- a correct answer to
        # a different question, not a missing input.
        from ...materials import available

        coating = QComboBox()
        coating.addItem("(none)", userData="")
        coating.setToolTip(
            "No coating gives relative efficiency -- a correct answer to a "
            "different question, not a missing input. Naming a material "
            "applies its Fresnel reflectivity and makes the result absolute; "
            "the provenance panel says which you got."
        )
        for name in available():
            coating.addItem(name, userData=name)
        coating.activated.connect(lambda _i: self.changed.emit())
        self._fields["coating"] = coating
        surface_form.addRow("Coating", coating)

        entry(surface_form, "Roughness σ (nm)", "roughness")

        scan = QGroupBox("Wavelengths (nm)")
        scan_form = QFormLayout(scan)
        entry(scan_form, "Start", "wavelength_start")
        entry(scan_form, "Stop", "wavelength_stop")
        entry(scan_form, "Points", "wavelength_count")

        for group in (grating, mount, surface, scan):
            column.addWidget(group)
        column.addStretch(1)

    def _check_fields_match_formstate(self) -> None:
        """Refuse to open rather than fail later on one code path.

        `read_form` builds a `FormState` by keyword from this dict, so a
        renamed or forgotten field is a TypeError at solve time -- visible
        only to whoever presses Solve. Checking here makes it immediate and
        says which field.
        """
        declared = {f.name for f in dataclasses.fields(FormState)}
        if self._fields.keys() != declared:
            raise AssertionError(
                "geometry fields do not match FormState: "
                f"missing {sorted(declared - self._fields.keys())}, "
                f"extra {sorted(self._fields.keys() - declared)}"
            )

    def _on_profile_change(self) -> None:
        """Show only the parameters the selected profile actually uses."""
        needed = PROFILE_FIELDS.get(
            self._fields["profile_kind"].currentText(), frozenset()
        )
        for key, widget in self._profile_rows.items():
            # setRowVisible hides the label with the field. Hiding the field
            # alone would leave a caption for a control that is not there.
            self._grating_form.setRowVisible(widget, key in needed)
        from_file = "profile_path" in needed
        self._load_button.setVisible(from_file)
        self._path_label.setVisible(from_file)

    def _on_mount_change(self) -> None:
        """Relabel the two angle fields; the mount decides what they mean."""
        primary, secondary = ANGLE_LABELS[self._fields["mount"].currentText()]
        self._angle_labels["alpha"].setText(primary)
        self._angle_labels["alpha"].show()
        self._fields["alpha"].show()

        if secondary is None:
            self._angle_labels["gamma"].hide()
            self._fields["gamma"].hide()
        else:
            self._angle_labels["gamma"].setText(secondary)
            self._angle_labels["gamma"].show()
            self._fields["gamma"].show()

    def read_form(self) -> FormState:
        return FormState(**{k: _value(w) for k, w in self._fields.items()})

    def load_profile(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load boundary profile", "", "PCGrate boundary (*.ggp);;All files (*)"
        )
        if path:
            self._fields["profile_path"].setText(path)
            self._path_label.setText(Path(path).name)
            self.changed.emit()
            self.solve_requested.emit()


def _value(widget: QWidget) -> str:
    """The *value* of an input, whichever kind it is.

    A combo's value is its `userData` when it has one, not its label. The
    coating combo shows "(none — relative efficiency)" for the empty choice,
    and `FormState` wants "" -- putting the label in would make the empty
    coating an unknown material name.
    """
    if isinstance(widget, QComboBox):
        data = widget.currentData()
        return widget.currentText() if data is None else str(data)
    return widget.text()
