# File formats

Three input formats are supported: the Bruker/Veeco Nanoscope `.spm` file,
the Asylum Research `.ibw` file, and a text export from Gwyddion. On the same
scan the `.spm` and text routes agree to within 0.002°. The `.ibw` route has
not yet been checked against a file from an instrument; see its section.

## Text exports

A text export is a tab-separated grid of heights in metres, preceded by a short
comment header:

```
# Channel: ZSensor
# Width: 2.000 µm
# Height: 1.566 µm
# Value units: m
```

The scan width is read from that header when it is present; otherwise the
fallback value in `SCAN_X_SIZE` is used.

Until recently this was the only supported input, and every file in `data/` is
a text export. Its drawback is that the export is an untracked manual step:
nothing in the file records which channel was exported, which scan direction was
used, or what processing was applied beforehand. The `_flatten` in those
filenames is the only indication that any processing was applied.

## Nanoscope `.spm`

This is the file the microscope writes. It consists of a latin-1 ASCII header of
`\Key: value` lines followed by the binary image planes.

One file holds several planes. A typical PeakForce scan contains four:

| Channel | Direction |
|---|---|
| Height Sensor | Retrace |
| Height Sensor | Trace |
| Peak Force Error | Retrace |
| Peak Force Error | Trace |

The header section for each plane gives the offset at which its bytes start
(`Data offset`), the width of each sample (`Bytes/pixel`, 2 or 4), and its
dimensions (`Samps/line` for columns, `Number of lines` for rows).

Rows and columns are read separately, because they are frequently unequal: the
bundled fixture is 512 wide by 401 lines, since that scan stopped early. Assuming
a square image would read past the end of the plane and reshape the next
channel's bytes into the height map.

### Which plane gets analysed

The plane is selected by `SPM_CHANNEL` and `SPM_DIRECTION`, which default to
Height Sensor / Retrace.

Only Height Sensor is offered in the GUI. Peak Force Error is a feedback error
signal and does not represent topography; passing it through a blaze-angle
pipeline produces a plausible-looking number with no physical meaning.

Retrace is the default because it is the plane from which this project's
existing Gwyddion exports were taken: it correlates 0.98 with the export, while
Trace correlates 0.90. With this default the two input routes agree by
construction. Both directions image the same grooves and either is a legitimate
choice, and the GUI offers a selector when a `.spm` is loaded.

### The height scale, and how to misread it

Raw samples are integers. They are converted to nanometres using a line such
as:

```
\@2:Z scale: V [Sens. ZsensSens] (0.0000000001872337 V/LSB) 0.804163 V
```

which gives `height = raw × (V/LSB) × Sens`.

The name in brackets must be used exactly as written. A file typically defines
both of the following:

```
\@Sens. Zsens:     V  32.46000 nm/V
\@Sens. ZsensSens: V 166.6319  nm/V
```

`Zsens` resembles `ZsensSens` with a redundant suffix removed, and using it
produces heights 5.13× too small. Nothing in the result appears wrong: the image
has the right shape, the grooves are in the right places, and the flattening
succeeds, but every reported blaze angle is incorrect.

This case is guarded by a test that compares the peak-to-peak height against the
known export and names the mistake in its failure message.

## Asylum `.ibw`

Asylum's MFP-3D and Cypher software is written in Igor Pro, and saves a scan as
an Igor binary wave (version 5): one 3-D wave whose layers are the channels,
named in the labels of its third dimension (`HeightRetrace`, `ZSensorTrace`,
`AmplitudeRetrace`, `PhaseRetrace`, ...), with the scan parameters in a text
note. The reader is `io/ibw.py`, written to WaveMetrics Technical Note PTN003,
and like the `.spm` reader it needs no library and recognises the file by its
content, in either byte order.

**Which channel.** The scan direction is the same setting as for `.spm`
(Retrace by default). The channel is chosen automatically: `ZSensor`, the
closed-loop sensor reading and the counterpart of Nanoscope's Height Sensor,
when the scan recorded one, and `Height` (from the piezo drive) otherwise. The
window shows which it used. Amplitude, phase and deflection are not topography
and are never chosen; `read_ibw(path, channel=...)` can name another channel
from Python.

**Scale and units.** Asylum stores heights in metres and sets the wave's X
scaling to metres per point, so heights are used as stored and the scan width
is the X scaling times the number of points, falling back to the note's
`ScanSize` when the wave has no physical scaling. A channel declared in other
units, or whose values reach a millimetre, is refused rather than read: that is
the `.ibw` form of the `.spm` scale trap, a constant factor on every height with
nothing visibly wrong.

**Not yet checked against an instrument file.** None was available when the
reader was written (2026-10-24). The tests build files to the specification,
including a synthetic blazed grating whose 25° blaze and 400 nm period come back
through the full analysis as 24.81° and 400.29 nm, which rules out a transposed
or rescaled read. They cannot rule out a misreading of the specification shared
by the reader and the tests. The first real Asylum scan, with an Igor or
Gwyddion export of the same scan, should be added as a fixture, as the `.spm`
reader's was.

## The two routes compared

The same scan was processed by both routes:

| Route | N | Mean blaze angle | σ | Period |
|---|---|---|---|---|
| Gwyddion `.txt` | 100 | 27.913° | 2.127° | 314.33 nm |
| Raw `.spm` | 100 | 27.911° | 2.126° | 314.33 nm |

The means differ by 0.002°, against a standard deviation of 2.13°.

The results are not bit-identical, and they are not expected to be: the export
was flattened by Gwyddion while the `.spm` is raw, and the two differ by 7 nm
RMS. This software, however, flattens every row-group profile itself before
fitting, which removes the same background. The Gwyddion flattening step is
therefore redundant for this pipeline.

A further finding was made during this verification: the `_flatten.0_00003`
companion file beside the raw scan is byte-identical to it. Gwyddion never wrote
modified data, and the flattening existed only in the text export.

## Practical notes

- The format is detected from the file's content and not from its extension.
  Nanoscope writes companion files with no extension at all, and these load
  correctly.
- Both routes return heights in metres, so no downstream step needs to know
  which format the data came from.
- A raw `.spm` occupies a few megabytes, roughly twice the size of the text
  export, and it carries every channel instead of one.
