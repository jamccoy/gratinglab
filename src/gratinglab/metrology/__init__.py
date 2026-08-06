"""
AFM Blaze Angle Analysis Package

A comprehensive toolkit for analyzing atomic force microscopy (AFM) data
of diffraction gratings to measure blaze angles.
"""

__version__ = '1.0.0'
__author__ = 'Jake McCoy'

# Import main workflow functions for easy access
from .workflows import (
    run_single_file_analysis,
    run_multiple_file_analysis,
    run_comparison_analysis
)

# Import analyzer for advanced usage
from .analyzer import analyze_single_file

__all__ = [
    'run_single_file_analysis',
    'run_multiple_file_analysis',
    'run_comparison_analysis',
    'analyze_single_file'
]