"""
Data management and aggregation.

This package combines scans and groups them by temperature.
"""

from .aggregation import (
    group_by_temperature,
    extract_temperatures_for_output,
    combine_scans
)

__all__ = [
    'group_by_temperature',
    'extract_temperatures_for_output',
    'combine_scans'
]