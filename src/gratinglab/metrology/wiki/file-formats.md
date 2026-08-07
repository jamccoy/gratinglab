# File formats

Two ways in: the instrument's own **Nanoscope `.spm`** file, or a **text export**
from Gwyddion. Both are supported, and on the same scan they agree to within
0.002°.

## Text exports

A tab-separated grid of heights in metres, with a short comment header:

```
# Channel: ZSensor
# Width: 2.000 µm
# Height: 1.566 µm
# Value units: m
```

The scan width is read from that header when present, otherwise the fallback in
`SCAN_X_SIZE` is used.

This was the only supported input until recently, and every file in `data/` is
one. The drawback is that the export is an untracked manual step: nothing in the
file records which channel was exported, which scan direction, or what processing
was applied first. The `_flatten` in those filenames is the only clue that
anything was.

## Nanoscope `.spm`

The file the microscope writes. A latin-1 ASCII header of `\Key: value` lines,
then the binary image planes.

One file holds **several planes**. A typical PeakForce scan contains four:

| Channel | Direction |
|---|---|
| Height Sensor | Retrace |
| Height Sensor | Trace |
| Peak Force Error | Retrace |
| Peak Force Error | Trace |

Each plane's header section says where its bytes start (`Data offset`), how wide
each sample is (`Bytes/pixel`, 2 or 4), and its dimensions (`Samps/line` for
columns, `Number of lines` for rows).

**Rows and columns are read separately.** They are frequently unequal — the
bundled fixture is 512 wide by 401 lines, because that scan stopped early.
Assuming a square image would read past the plane and reshape the next channel's
bytes into the height map.

### Which plane gets analysed

`SPM_CHANNEL` and `SPM_DIRECTION`, defaulting to **Height Sensor / Retrace**.

Only Height Sensor is offered in the GUI. Peak Force Error is a feedback error
signal, not topography; running it through a blaze-angle pipeline produces a
confident-looking number that means nothing.

Retrace is the default because it is the plane this project's existing Gwyddion
exports were taken from — it correlates 0.98 with the export, where Trace
correlates 0.90. That makes the two input routes agree by construction rather
than by luck. Both directions image the same grooves and either is a legitimate
choice; the GUI offers a selector when a `.spm` is loaded.

### The height scale, and the trap in it

Raw samples are integers. Converting them to nanometres uses a line like:

```
\@2:Z scale: V [Sens. ZsensSens] (0.0000000001872337 V/LSB) 0.804163 V
```

giving `height = raw × (V/LSB) × Sens`.

**The name in brackets must be used verbatim.** A file typically defines both:

```
\@Sens. Zsens:     V  32.46000 nm/V
\@Sens. ZsensSens: V 166.6319  nm/V
```

`Zsens` looks like `ZsensSens` with a redundant suffix stripped, and reaching for
it produces heights **5.13× too small**. Nothing looks wrong: the image is the
right shape, the grooves sit in the right places, the flattening succeeds, and
every reported blaze angle is quietly incorrect.

This is guarded by a test that compares peak-to-peak against the known export and
names the trap in its failure message.

## The two routes compared

Same scan, both ways:

| Route | N | Mean blaze angle | σ | Period |
|---|---|---|---|---|
| Gwyddion `.txt` | 100 | 27.913° | 2.127° | 314.33 nm |
| Raw `.spm` | 100 | 27.911° | 2.126° | 314.33 nm |

A difference of **0.002°**, against a standard deviation of 2.13°.

They are not bit-identical, and should not be: the export was flattened by
Gwyddion, the `.spm` is raw, and the two differ by 7 nm RMS. But this software
flattens every row-group profile itself before fitting, which removes the same
background. The Gwyddion step turns out to be redundant for this pipeline.

One incidental finding while verifying that: the `_flatten.0_00003` companion
file sitting beside the raw scan is **byte-identical** to it. Gwyddion never
wrote modified data — the flattening existed only in the text export.

## Practical notes

- Format is detected by **content**, not extension. Nanoscope writes companion
  files with no extension at all, and those load fine.
- Heights are returned in metres from both routes, so nothing downstream needs to
  know which format it came from.
- A raw `.spm` is a few megabytes against roughly half that for the text export,
  and carries every channel rather than one.
