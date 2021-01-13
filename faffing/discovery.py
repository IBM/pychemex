import os
import pandas as pd

data_dir = "../data"
raw_file = "/patents_ochem_enamine_bradley_begstrom_training_.csv"
os_agnostic_file_location = os.path.normpath(f"{data_dir}{raw_file}")

raw_data = pd.read_csv(
    os_agnostic_file_location,
    error_bad_lines=False,
    usecols=["SMILES", "NAME"]
)

raw_data.iloc[:10, :].to_csv(
    os.path.normpath(f"{data_dir}/smiles_only_ten.smi"),
    sep=" ",
    index=False,
    header=False
)

test = None
