"""File formats other programs read or write.

Each module handles one format, and nothing is imported here, so a caller
names the module it wants (``from gratinglab.io.ggp import read_ggp``):

``ggp``
    PCGrate's ``.ggp`` groove profiles, read in all five border types and
    written in the polygonal one.
``boundary``
    A measured groove profile saved with its period and provenance, as
    gratinglab's own ``.json``, which reads back, or a write-only ``.csv``.
    A ``.ggp`` cannot carry the period, which is why this exists beside it.
``efficiency_table``
    Efficiency tables exported by PCGrate-SX, read for comparison.
``txl``
    Raith ``.txl`` layouts of a grating, written for electron-beam writing.
"""
