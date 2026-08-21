"""
Metrology test setup.

The backend and the Qt environment are settled once, in `tests/conftest.py`, so
that the two suites cannot disagree about them. What remains here is the `src/`
path insert, kept only so these files still run directly as scripts -- several
carry an `if __name__ == '__main__'` runner. Under pytest it is redundant with
`pythonpath = ["src"]` in pyproject.toml, and harmless.
"""
import os
import sys

SRC = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    'src')
if SRC not in sys.path:
    sys.path.insert(0, SRC)
