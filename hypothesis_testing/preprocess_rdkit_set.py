import pandas as pd
from typing import List, Tuple


def load(data_path: str) -> Tuple[pd.DataFrame, List, List]:
    data = pd.read_csv(data_path)
    # remove some weird stuff
    data = data[data["qed"] != -666]

    # get fragment columns
    frag_cols = []
    for col in list(data.columns):
        if col[:3] == "fr_":
            frag_cols.append(col)

    # get features that are not fragments
    features = []
    for col in data.columns:
        if col not in frag_cols and col not in ["SMILES", "MP"]:
            features.append(col)

    return data, frag_cols, features
