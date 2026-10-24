"""The grating geometry: period, profile, mount, wavelengths.

This panel was extracted from the original `_build_inputs` of
`main_window.py`, without scalar's own "Scalar solver" groupbox and its
Solve/Cancel/Export controls. Those are specific to a solver and remain in
that solver's own tab (`ScalarTab`, and eventually whatever tab RCWA
requires). This panel is shared and built once, not once per solver, because
a `Problem` and an `Illumination` mean the same thing whichever solver is
about to consume them.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QPalette
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
    """Grating, mount and wavelength inputs, shown in the persistent left sidebar.

    The panel is not wrapped in its own scroll area, because the way it is
    embedded (at present, stacked above a solver's own controls under one
    shared scroller) is for the caller to decide.
    """

    #: Pressing Enter in a field, or finishing a "Load profile…" pick. This
    #: panel does not know what a solve is, only that something changed
    #: enough that one might be wanted; `MainWindow` connects it to `solve()`.
    solve_requested = Signal()

    #: The form was edited *by a user*. Distinct from
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
        # The palette's muted role rather than a fixed grey, so it follows
        # the theme (a stylesheet colour would also freeze the label's palette).
        self._path_label.setForegroundRole(QPalette.ColorRole.Mid)
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
        coating = QComboBox()
        coating.setToolTip(
            "No coating gives relative efficiency -- a correct answer to a "
            "different question, not a missing input. Naming a material "
            "makes the result absolute: the Scalar tab applies its Fresnel "
            "reflectivity, and the Integral tab reads its complex index when "
            'Conductivity is "tabulated" (and ignores it when "perfect"). '
            "The provenance panel says which you got."
        )
        coating.activated.connect(lambda _i: self.changed.emit())
        self._fields["coating"] = coating
        self.refresh_coatings()
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
        self._forms = (grating_form, mount_form, surface_form, scan_form)
        self._align_forms()

    def _align_forms(self) -> None:
        """Give the four groups one label column, so the inputs line up.

        Each group lays out its own form, so each sized its label column to
        its own longest label and the input boxes started at four different
        offsets down the dock (64 px apart at the widest). The shared width
        is measured over every text a label can take, the mount labels'
        alternatives included, so switching mount never moves the inputs.
        """
        metrics = self.fontMetrics()
        texts = [
            text
            for pair in ANGLE_LABELS.values()
            for text in pair
            if text is not None
        ]
        labels: list[QLabel] = []
        for form in self._forms:
            form.setFieldGrowthPolicy(
                QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
            )
            form.setLabelAlignment(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
            )
            for row in range(form.rowCount()):
                item = form.itemAt(row, QFormLayout.ItemRole.LabelRole)
                label = item.widget() if item is not None else None
                if isinstance(label, QLabel):
                    labels.append(label)
                    texts.append(label.text())
        # A little slack over the text, as a QLabel's own size hint carries.
        width = max(metrics.horizontalAdvance(t) for t in texts) + 6
        for label in labels:
            label.setFixedWidth(width)

    def refresh_coatings(self) -> None:
        """Re-read the material library into the Coating dropdown.

        The choice is kept when its material is still installed. When it is
        not, the dropdown falls back to "(none)" and emits ``changed``, since
        the problem the form describes has changed under the user. A library
        that cannot be read lists no materials here; the Setup tab says why.
        """
        from ...materials import available

        combo = self._fields["coating"]
        current = combo.currentData() or ""
        try:
            names = available()
        except ValueError:  # two tables for one name; Setup reports it
            names = ()
        combo.blockSignals(True)
        try:
            combo.clear()
            combo.addItem("(none)", userData="")
            for name in names:
                combo.addItem(name, userData=name)
            index = combo.findData(current)
            combo.setCurrentIndex(max(index, 0))
        finally:
            combo.blockSignals(False)
        if current and index < 0:
            self.changed.emit()

    def _check_fields_match_formstate(self) -> None:
        """Refuse to open rather than fail later on one code path.

        `read_form` builds a `FormState` by keyword from this dict, so a
        renamed or forgotten field would otherwise cause a TypeError at solve
        time, visible only to whoever presses Solve. Checking here reports the
        problem immediately and names the field.
        """
        declared = {f.name for f in dataclasses.fields(FormState)}
        if self._fields.keys() != declared:
            raise AssertionError(
                "geometry fields do not match FormState: "
                f"missing {sorted(declared - self._fields.keys())}, "
                f"extra {sorted(self._fields.keys() - declared)}"
            )

    def _on_profile_change(self) -> None:
        """Show only the parameters the selected profile uses."""
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
            self, "Load boundary profile", "",
            "Boundary profiles (*.json *.ggp);;"
            "gratinglab boundary, with period (*.json);;"
            "PCGrate boundary (*.ggp);;All files (*)",
        )
        if path:
            self.load_boundary_file(path)

    def load_boundary_file(self, path: "str | Path") -> None:
        """Use a boundary file as the profile, and its period if it has one.

        A gratinglab ``.json`` carries the period it was measured at, so the
        period field is set from it; a ``.ggp`` cannot, so the field is left as
        it is, which is the one case where the user still has to supply it.
        The metrology window's "Send to efficiency window" also calls this
        method.
        """
        path = Path(path)
        note = ""
        if path.suffix.lower() == ".json":
            from ...io.boundary import read_boundary_json

            try:
                saved = read_boundary_json(path)
            except (OSError, ValueError):
                pass  # the form's own validation reports the reason
            else:
                self._fields["period"].setText(f"{saved.period_nm:.6g}")
                note = f" (period {saved.period_nm:.6g} nm from the file)"
        kind = self._fields["profile_kind"]
        kind.setCurrentIndex(kind.findText("From file"))
        self._fields["profile_path"].setText(str(path))
        self._path_label.setText(path.name + note)
        self.changed.emit()
        self.solve_requested.emit()


def _value(widget: QWidget) -> str:
    """The value of an input, whichever kind it is.

    A combo's value is its `userData` when it has one, not its label. The
    coating combo shows "(none)" for the empty choice, while `FormState`
    expects "", and passing the label would turn the empty coating into an
    unknown material name.
    """
    if isinstance(widget, QComboBox):
        data = widget.currentData()
        return widget.currentText() if data is None else str(data)
    return widget.text()
