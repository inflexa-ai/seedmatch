"""The dtype of the text columns in the tables that seedmatch returns."""

from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd
from pandas.api.extensions import ExtensionDtype


def text_dtype() -> np.dtype[Any] | ExtensionDtype:
    """Return the dtype that the installed pandas gives a column of text.

    The dtype is `object` on pandas 2, and `str` on pandas 3 or with the option
    `future.infer_string`.
    """
    return pd.Series([""]).dtype


def with_text_dtype(frame: pd.DataFrame, names: Iterable[str]) -> pd.DataFrame:
    """Return the frame with the dtype of `text_dtype` for each column of `names` of only text.

    pandas 3 gives `str` to a column of text, but `object` to an empty column or to a column of
    only missing values. Thus the dtype of a table can change with its rows. A column that has the
    dtype already does not change. A column with a value that is not text and not missing, for
    example an int ID of the caller, keeps its dtype and its values.
    """
    dtype = text_dtype()
    changes = {
        name: dtype
        for name in names
        if frame[name].dtype != dtype
        and pd.api.types.infer_dtype(frame[name], skipna=True) in {"string", "empty"}
    }
    return frame.astype(changes) if changes else frame
