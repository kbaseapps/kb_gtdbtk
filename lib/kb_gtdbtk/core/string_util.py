"""
Contains utility functions for creating and manipulating strings for use in multiple modules.
"""
from datetime import datetime

def now_ISOish() -> str:
    """
    Returns the current time in a format roughly resembling ISO:
    YYYYMMDD_HHMMSS
    I.e.: 20151112_162213
    """
    return datetime.now().strftime("%Y%m%d_%H%M%S")
