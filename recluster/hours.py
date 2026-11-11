"""Clock-time parsing shared by the claim spec and the time-of-day check."""

from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd


def parse_hours(values: pd.Series) -> pd.Series:
    """Clock times as hours in [0, 24): accepts numbers, 'HH:MM', 'HH:MM:SS' and datetimes."""
    out = []
    for v in values:
        if v is None or (isinstance(v, float) and math.isnan(v)):
            out.append(np.nan)
            continue
        if isinstance(v, (int, float, np.integer, np.floating)):
            out.append(float(v) % 24.0)
            continue
        if isinstance(v, (pd.Timestamp, np.datetime64)):
            ts = pd.Timestamp(v)
            out.append(ts.hour + ts.minute / 60 + ts.second / 3600)
            continue
        s = str(v).strip()
        m = re.fullmatch(r"(\d{1,2}):(\d{2})(?::(\d{2}))?", s)
        if m:
            h, mi, se = int(m.group(1)), int(m.group(2)), int(m.group(3) or 0)
            out.append((h + mi / 60 + se / 3600) % 24.0)
            continue
        try:
            out.append(float(s) % 24.0)
        except ValueError:
            try:
                ts = pd.Timestamp(s)
                out.append(ts.hour + ts.minute / 60 + ts.second / 3600)
            except (ValueError, TypeError):
                out.append(np.nan)
    return pd.Series(out, index=values.index, dtype=float)
