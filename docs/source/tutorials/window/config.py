"""Configuration module for the window convolution tutorials.

This module defines the root directory and data directory for the tutorials.
"""
from pathlib import Path

ROOTDIR = Path(__file__).parent.parent
DATADIR = ROOTDIR / "data"
