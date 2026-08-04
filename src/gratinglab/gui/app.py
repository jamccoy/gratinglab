"""The window.

Deliberately thin. Every number shown here comes from
:func:`gratinglab.compare.sweep`; every validation comes from
:mod:`gratinglab.gui.state`. Nothing in this module computes physics, so a bug
here misplaces a widget rather than producing a wrong answer.

Run with ``gratinglab-gui``.
"""

from __future__ import annotations

import base64
import csv
import re
import sys
from pathlib import Path

from . import mathtext, provenance
from .docs import general_pages, theory_pages
from .state import (
    ANGLE_LABELS,
    MOUNTS,
    PROFILE_FIELDS,
    PROFILE_KINDS,
    FormErrors,
    FormState,
    build,
)

#: `## Heading` / `### Sub-heading` -- the only markdown structure the theory
#: viewer styles beyond math. Tables and bullet lists stay literal; see
#: mathtext.py's module docstring for why that scope is deliberate.
_HEADING_RE = re.compile(r"^(#{2,6})\s+(.*)$")
_BOLD_RE = re.compile(r"\*\*(.+?)\*\*")

#: Display equations render a little larger than inline ones, roughly matching
#: how a printed page distinguishes a headline equation from a symbol in running
#: text.
_DISPLAY_DPI = 170
_INLINE_DPI = 130

_TK_HELP = """\
gratinglab-gui needs Tk, which this Python was built without.

  macOS/Homebrew:  brew install python-tk@3.12
  pyenv:           install a Python built with tcl-tk
  Debian/Ubuntu:   sudo apt install python3-tk

Everything except the GUI works without it.
"""


def _require_tk():
    """Import tkinter with an explanation instead of a traceback."""
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except ImportError as exc:  # pragma: no cover - depends on the build
        raise SystemExit(_TK_HELP) from exc
    return tk, ttk, filedialog, messagebox


class GratingLabApp:
    """Single-window explorer for one solver over one wavelength scan."""

    def __init__(self, root) -> None:
        tk, ttk, filedialog, messagebox = _require_tk()
        self._tk, self._ttk = tk, ttk
        self._filedialog, self._messagebox = filedialog, messagebox

        self.root = root
        self.root.title("GratingLab")
        self.root.geometry("1180x780")

        self._vars: dict[str, object] = {}
        self._scan = None

        self._build_menu()

        body = ttk.Frame(root, padding=8)
        body.pack(fill="both", expand=True)

        self._build_inputs(ttk.Frame(body, padding=(0, 0, 12, 0)))
        self._build_plots(body)
        self._build_provenance(root)

        self._on_mount_change()
        self._on_profile_change()
        self.solve()

    # -- layout ----------------------------------------------------------

    def _build_inputs(self, frame) -> None:
        tk, ttk = self._tk, self._ttk
        frame.pack(side="left", fill="y")
        row = 0

        def heading(text):
            nonlocal row
            ttk.Label(frame, text=text, font=("TkDefaultFont", 11, "bold")).grid(
                row=row, column=0, columnspan=2, sticky="w", pady=(10, 2)
            )
            row += 1

        def entry(label, key, default):
            nonlocal row
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w")
            var = tk.StringVar(value=default)
            self._vars[key] = var
            ttk.Entry(frame, textvariable=var, width=12).grid(
                row=row, column=1, sticky="ew", pady=1
            )
            row += 1
            return var

        def combo(label, key, values, default, command=None):
            nonlocal row
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w")
            var = tk.StringVar(value=default)
            self._vars[key] = var
            box = ttk.Combobox(
                frame, textvariable=var, values=list(values), width=10, state="readonly"
            )
            box.grid(row=row, column=1, sticky="ew", pady=1)
            if command:
                box.bind("<<ComboboxSelected>>", lambda _e: command())
            row += 1
            return var

        defaults = FormState()

        heading("Grating")
        entry("Period (nm)", "period", defaults.period)
        combo("Profile", "profile_kind", PROFILE_KINDS, defaults.profile_kind,
              self._on_profile_change)

        self._profile_rows: dict[str, list] = {}
        for key, label, default in (
            ("blaze_angle", "Blaze δ (deg)", defaults.blaze_angle),
            ("antiblaze_angle", "Anti-blaze (deg)", defaults.antiblaze_angle),
            ("depth_fraction", "Depth / period", defaults.depth_fraction),
            ("duty_cycle", "Duty cycle", defaults.duty_cycle),
        ):
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w")
            var = tk.StringVar(value=default)
            self._vars[key] = var
            widget = ttk.Entry(frame, textvariable=var, width=12)
            widget.grid(row=row, column=1, sticky="ew", pady=1)
            self._profile_rows[key] = [frame.grid_slaves(row=row, column=0)[0], widget]
            row += 1

        self._vars["profile_path"] = tk.StringVar(value="")
        self._path_label = ttk.Label(frame, text="", foreground="#555", wraplength=170)
        self._path_label.grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1
        self._load_button = ttk.Button(
            frame, text="Load profile…", command=self.load_profile
        )
        self._load_button.grid(row=row, column=0, columnspan=2, sticky="ew", pady=2)
        row += 1

        heading("Mount")
        combo("Geometry", "mount", MOUNTS, defaults.mount, self._on_mount_change)
        self._angle_labels = {}
        for key, default in (("alpha", defaults.alpha), ("gamma", defaults.gamma)):
            label = ttk.Label(frame, text="")
            label.grid(row=row, column=0, sticky="w")
            var = tk.StringVar(value=default)
            self._vars[key] = var
            widget = ttk.Entry(frame, textvariable=var, width=12)
            widget.grid(row=row, column=1, sticky="ew", pady=1)
            self._angle_labels[key] = (label, widget)
            row += 1

        heading("Wavelengths (nm)")
        entry("Start", "wavelength_start", defaults.wavelength_start)
        entry("Stop", "wavelength_stop", defaults.wavelength_stop)
        entry("Points", "wavelength_count", defaults.wavelength_count)

        heading("Scalar solver")
        entry("Quadrature pts", "quadrature_points", defaults.quadrature_points)

        ttk.Separator(frame).grid(row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1
        ttk.Button(frame, text="Solve", command=self.solve).grid(
            row=row, column=0, columnspan=2, sticky="ew"
        )
        row += 1
        ttk.Button(frame, text="Export CSV…", command=self.export_csv).grid(
            row=row, column=0, columnspan=2, sticky="ew", pady=2
        )
        frame.columnconfigure(1, weight=1)

    def _build_plots(self, body) -> None:
        import matplotlib

        matplotlib.use("TkAgg")
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 6), dpi=100, layout="constrained")
        self._profile_axes = self._figure.add_subplot(2, 1, 1)
        self._efficiency_axes = self._figure.add_subplot(2, 1, 2)

        self._canvas = FigureCanvasTkAgg(self._figure, master=body)
        self._canvas.get_tk_widget().pack(side="right", fill="both", expand=True)

    def _build_provenance(self, root) -> None:
        tk, ttk = self._tk, self._ttk
        frame = ttk.Frame(root, padding=(8, 0, 8, 8))
        frame.pack(fill="x", side="bottom")
        ttk.Label(
            frame, text="Provenance", font=("TkDefaultFont", 11, "bold")
        ).pack(anchor="w")
        self._provenance = tk.Text(frame, height=7, wrap="word", relief="solid",
                                   borderwidth=1)
        self._provenance.pack(fill="x")
        for tag, color in provenance.TAG_COLORS.items():
            self._provenance.tag_configure(tag, foreground=color)
        self._provenance.configure(state="disabled")

    def _build_menu(self) -> None:
        """Help menu: what each solver actually computes, and About.

        Deliberately the only menu. This is not a general-purpose application
        with File/Edit/View concerns -- it is one window with one job, and the
        single thing worth a menu for is explaining the math, which the window
        itself never does.
        """
        tk = self._tk
        menubar = tk.Menu(self.root)

        help_menu = tk.Menu(menubar, tearoff=False)
        help_menu.add_command(label="About GratingLab", command=self.show_about)
        help_menu.add_separator()

        # Foundational, solver-independent reference first (the generalized
        # grating equation, angle conventions), then one page per method.
        for page in general_pages():
            label = page.title + ("" if page.available else " (not written yet)")
            help_menu.add_command(
                label=label, command=lambda p=page: self.show_theory(p)
            )
        help_menu.add_separator()
        for page in theory_pages():
            label = f"{page.title} Theory" + ("" if page.available else " (not written yet)")
            help_menu.add_command(
                label=label, command=lambda p=page: self.show_theory(p)
            )
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.config(menu=menubar)

    # -- reactions -------------------------------------------------------

    def _on_profile_change(self) -> None:
        """Show only the parameters the selected profile actually uses."""
        kind = self._vars["profile_kind"].get()
        needed = PROFILE_FIELDS.get(kind, frozenset())
        for key, widgets in self._profile_rows.items():
            for widget in widgets:
                widget.grid() if key in needed else widget.grid_remove()
        from_file = kind == "From file"
        self._load_button.grid() if from_file else self._load_button.grid_remove()
        self._path_label.grid() if from_file else self._path_label.grid_remove()

    def _on_mount_change(self) -> None:
        """Relabel the two angle fields; the mount decides what they mean."""
        primary, secondary = ANGLE_LABELS[self._vars["mount"].get()]
        label, widget = self._angle_labels["alpha"]
        label.configure(text=primary)
        label.grid(), widget.grid()
        label, widget = self._angle_labels["gamma"]
        if secondary is None:
            label.grid_remove(), widget.grid_remove()
        else:
            label.configure(text=secondary)
            label.grid(), widget.grid()

    def _read_form(self) -> FormState:
        return FormState(
            **{key: var.get() for key, var in self._vars.items()}  # type: ignore[arg-type]
        )

    # -- actions ---------------------------------------------------------

    def load_profile(self) -> None:
        path = self._filedialog.askopenfilename(
            title="Load boundary profile",
            filetypes=[("PCGrate boundary", "*.ggp"), ("All files", "*")],
        )
        if path:
            self._vars["profile_path"].set(path)
            self._path_label.configure(text=Path(path).name)
            self.solve()

    def solve(self) -> None:
        """Validate, solve, plot. All physics happens in the core."""
        from ..checks import check_energy_balance
        from ..compare import sweep

        try:
            parsed = build(self._read_form())
        except FormErrors as exc:
            self._show_errors(exc)
            return

        scan = sweep(
            parsed.problem,
            parsed.illumination,
            parsed.wavelengths,
            ["scalar"],
            options={"scalar": parsed.options},
        )[0]
        self._scan = scan

        self._draw_profile(parsed)
        self._draw_efficiency(scan)
        self._canvas.draw_idle()
        self._show_provenance(scan, check_energy_balance(scan), parsed)

    def export_csv(self) -> None:
        if self._scan is None:
            self._messagebox.showinfo("Nothing to export", "Solve first.")
            return
        path = self._filedialog.asksaveasfilename(
            title="Export efficiencies", defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
        )
        if not path:
            return
        rows = self._scan.to_records()
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        self._messagebox.showinfo("Exported", f"{len(rows)} rows to {Path(path).name}")

    def show_about(self) -> None:
        from .. import __version__

        self._messagebox.showinfo(
            "About GratingLab", provenance.about_text(__version__)
        )

    def show_theory(self, page) -> None:
        """Open the read-only viewer for one theory or reference page.

        Math is typeset -- rasterised by matplotlib's mathtext and embedded as
        images -- rather than shown as literal `$...$` source. Headings and
        `**bold**` get light styling; tables and bullet lists stay literal
        text, since real column/list layout in a `Text` widget would be a lot
        of engineering for content that is already readable as-is (see
        mathtext.py's module docstring). A non-rigorous solver gets a banner
        up top, matching the practice of surfacing an approximation rather
        than hiding it.
        """
        tk = self._tk
        window = tk.Toplevel(self.root)
        window.title(page.title)
        window.geometry("820x680")
        # PhotoImage instances are garbage-collected as soon as nothing in
        # Python holds a reference, even though Tk is still displaying them --
        # this list is that reference, and lives exactly as long as the window.
        window._photo_refs = []  # type: ignore[attr-defined]

        text = tk.Text(window, wrap="word", font=("Menlo", 12), padx=14, pady=12)
        scrollbar = tk.Scrollbar(window, command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        # Scrollbar packed first: it claims its strip before Text expands to
        # fill whatever remains, rather than being squeezed out afterwards.
        scrollbar.pack(side="right", fill="y")
        text.pack(side="left", fill="both", expand=True)

        text.tag_configure("banner", foreground="#a5370d", font=("Menlo", 12, "bold"))
        text.tag_configure("h2", font=("Menlo", 15, "bold"), spacing1=8, spacing3=6)
        text.tag_configure("h3", font=("Menlo", 13, "bold"), spacing1=6, spacing3=4)
        text.tag_configure("bold", font=("Menlo", 12, "bold"))
        text.tag_configure("center", justify="center")
        text.tag_configure("unrendered", foreground="#a5370d", font=("Menlo", 11))

        if not page.rigorous:
            text.insert(
                "end",
                "⚠ Approximate method — see Limits below.\n\n",
                "banner",
            )

        for segment in mathtext.split_segments(page.text):
            if segment.kind == "text":
                self._insert_text_segment(text, segment.content)
            else:
                self._insert_math_segment(text, window, segment)

        text.configure(state="disabled")

    def _insert_text_segment(self, text, content: str) -> None:
        """Insert a text segment, styling `## headings` and `**bold**`.

        Split on the segment's own newlines and rejoined by re-inserting `\n`
        only *between* lines -- never before the first or after the last --
        so a segment that starts or ends mid-line (because a math span sits
        right next to it) does not gain a newline that was not in the source.
        Headings are safe to detect per-line this way because, in every page
        this viewer shows today, a heading is always a complete line with no
        adjacent math span on it.
        """
        lines = content.split("\n")
        for index, line in enumerate(lines):
            heading = _HEADING_RE.match(line)
            if heading:
                level = len(heading.group(1))
                text.insert("end", heading.group(2), "h2" if level <= 2 else "h3")
            else:
                self._insert_bold_spans(text, line)
            if index != len(lines) - 1:
                text.insert("end", "\n")

    def _insert_bold_spans(self, text, line: str) -> None:
        position = 0
        for match in _BOLD_RE.finditer(line):
            text.insert("end", line[position : match.start()])
            text.insert("end", match.group(1), "bold")
            position = match.end()
        text.insert("end", line[position:])

    def _insert_math_segment(self, text, window, segment) -> None:
        """Render one math span to an image and embed it, or fall back to
        visibly-tagged raw source if mathtext could not parse it."""
        dpi = _DISPLAY_DPI if segment.display else _INLINE_DPI
        png = mathtext.render_math_png(segment.content, dpi=dpi)

        if png is None:
            delimiter = "$$" if segment.display else "$"
            text.insert(
                "end",
                f"{delimiter}{segment.content}{delimiter}",
                "unrendered",
            )
            return

        photo = self._tk.PhotoImage(
            data=base64.b64encode(png).decode("ascii"), format="png"
        )
        window._photo_refs.append(photo)  # type: ignore[attr-defined]

        if segment.display:
            text.insert("end", "\n")
            start = text.index("end-1c")
            text.image_create("end", image=photo)
            text.insert("end", "\n")
            text.tag_add("center", start, text.index("end-1c"))
        else:
            text.image_create("end", image=photo)

    # -- drawing ---------------------------------------------------------

    def _draw_profile(self, parsed) -> None:
        import numpy as np

        axes = self._profile_axes
        axes.clear()
        t = np.linspace(0.0, 1.0, 600, endpoint=False)
        # Drawn from the same Profile the solver integrates, so what is shown
        # is literally what is computed.
        height = parsed.problem.height_nm(np.concatenate([t, t + 1.0]))
        axes.plot(np.concatenate([t, t + 1.0]), height, color="#1f3b73", lw=1.8)
        axes.fill_between(
            np.concatenate([t, t + 1.0]), 0, height, color="#1f3b73", alpha=0.12
        )
        axes.set_xlabel("position / period  (two periods shown)")
        axes.set_ylabel("height (nm)")
        axes.set_title(
            f"{type(parsed.problem.profile).__name__} — "
            f"period {parsed.problem.period:g} nm, depth {parsed.problem.depth:.4g} nm",
            fontsize=10,
        )
        axes.grid(alpha=0.25)

    def _draw_efficiency(self, scan) -> None:
        axes = self._efficiency_axes
        axes.clear()
        for index, order in enumerate(scan.orders):
            values = scan.efficiency[:, index]
            if values.max() < 1e-4:
                continue
            axes.plot(scan.wavelengths, values, lw=1.3, label=f"m={int(order):+d}")
        axes.plot(scan.wavelengths, scan.total, "k--", lw=1.0, alpha=0.7, label="Σ")
        axes.axhline(1.0, color="#b00020", lw=0.8, ls=":", alpha=0.7)
        axes.set_xlabel("wavelength (nm)")
        axes.set_ylabel("efficiency")
        axes.set_ylim(bottom=0)
        axes.grid(alpha=0.25)
        handles, _ = axes.get_legend_handles_labels()
        if handles:
            axes.legend(fontsize=7, ncol=max(1, len(handles) // 8), loc="upper right")

    # -- provenance ------------------------------------------------------

    def _paint(self, lines) -> None:
        """Paint `provenance.Line` tuples into the panel.

        The whole of this layer's involvement. Which line carries which tag is
        decided in `gui/provenance.py`, where it can be tested without a
        window -- and where the rule that a correct default must not be styled
        as a problem is written down.
        """
        self._provenance.configure(state="normal")
        self._provenance.delete("1.0", "end")
        for line in lines:
            self._provenance.insert("end", line.text, line.tag or "")
        self._provenance.configure(state="disabled")

    def _show_errors(self, exc: FormErrors) -> None:
        self._paint(provenance.error_lines(exc.errors))

    def _show_provenance(self, scan, energy, parsed) -> None:
        self._paint(
            provenance.provenance_lines(scan, energy, parsed.lambda_over_period)
        )


def _bring_to_front(root) -> None:
    """Put the window in front on launch.

    A Tk window started from a shell stub opens *behind* whatever is frontmost
    on macOS, because the process has no bundle identity of its own. Raising it
    and briefly asserting topmost is the standard remedy; the attribute is
    released immediately so the window behaves normally afterwards.
    """
    root.lift()
    root.attributes("-topmost", True)
    root.after(200, lambda: root.attributes("-topmost", False))
    try:
        # Ask macOS to activate this process. Harmless if it fails.
        from subprocess import run

        run(
            ["osascript", "-e",
             'tell application "System Events" to set frontmost of the first '
             'process whose unix id is {} to true'.format(__import__("os").getpid())],
            capture_output=True, timeout=3,
        )
    except Exception:  # pragma: no cover - best effort only
        pass


def main() -> int:
    """Console-script entry point."""
    tk, _, _, _ = _require_tk()
    root = tk.Tk()
    GratingLabApp(root)
    _bring_to_front(root)
    root.mainloop()
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
