from __future__ import annotations

import uuid
from dataclasses import dataclass
import pandas as pd


@dataclass
class Precedent:
    name: str
    rxn_smi: str
    precedent_id: str
    data: dict
    similarity: float = 0


def make_precedents(df: pd.DataFrame, name_col: str, id_col: str):
    """ Take a dataframe of precedents and create precedent dataclass per row """
    precedents = []
    for i, row in df.iterrows():
        row_no_nan = row.dropna()
        precedents.append(make_single_precedent(row_no_nan, name_col, id_col))
    return precedents


def parse_row_to_dict(row: pd.Series) -> dict:
    data = dict(row)
    for key, value in data.items():
        # check if value is a str, float, or int
        if isinstance(value, str):
            continue
        elif isinstance(value, float):
            continue
        elif isinstance(value, int):
            continue
        else:  # otherwise make it a str
            data[key] = str(value)
    return data

def make_single_precedent(row: pd.Series | dict, name_col: str, id_col: str):
    data = parse_row_to_dict(row)

    precedent_id = data.get(id_col, uuid.uuid4())

    if data.get(name_col, None) is None:
        name = precedent_id
    else:
        name = row[name_col]

    if isinstance(name, list):
        name = ', '.join(name)

    return Precedent(name=name,
                     rxn_smi=row['rxn_smi'],
                     precedent_id=precedent_id,
                     data=data,
                     similarity=float(row['similarity']))
