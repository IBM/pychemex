from typing import AnyStr, Any, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors


class Molecule:
    """
    Stores SMILE and descriptor for a given structure. Currently only taking RDKit descriptors but possible to expand.
    """

    def __init__(self, smile: AnyStr, melting_point: float):
        """
        :param smile: string encoding the molecule information
        :param melting_point: the recorded melting_point in Kelvin
        """
        self.smile = smile
        self.melting_point = melting_point
        self.descriptors = dict()

    def add_descriptor(self, descriptor_name: AnyStr, descriptor: Any):
        """Add descriptor value and name"""
        self.descriptors[descriptor_name] = descriptor

    def add_descriptors_from_lists(self, descriptor_name: List[str], descriptor_value: List[Any]):
        """Adds descriptors from lists"""
        self.descriptors.update(zip(descriptor_name, descriptor_value))

    def add_basic_info_to_descriptors(self):
        """ Add the smiles and melting points to dictionary of descriptors"""
        self.descriptors['SMILES'] = self.smile
        self.descriptors['MP'] = self.melting_point


def read_mp_data(file_path: AnyStr) -> pd.DataFrame:
    """Reads melting point data"""
    return pd.read_csv(file_path)


def create_descriptor_df(molecules: List[Molecule]) -> pd.DataFrame:
    """Geneates dataframe from list of molecules"""
    return pd.DataFrame([mol.descriptors for mol in molecules])


def generate_rdkit_mol(SMILE: AnyStr) -> Chem.rdchem.Mol:
    """generates Rdkit molecule to carry out calculations with"""
    return Chem.MolFromSmiles(SMILE)


data = read_mp_data('D:\pychemex\data\melting_point\MP_dataset.csv')
molecules: List[Molecule] = []
list_of_descriptors = [desc[0] for desc in Descriptors.descList]
all_desc = MoleculeDescriptors.MolecularDescriptorCalculator(list_of_descriptors)
for index, row in data.iterrows():
    print(f"{index} out of {len(data)}")
    mol1 = Molecule(smile=row['SMILES'], melting_point=row['MP'])
    mol1_desc = all_desc.CalcDescriptors(generate_rdkit_mol(mol1.smile))
    mol1.add_descriptors_from_lists(list_of_descriptors, mol1_desc)
    mol1.add_basic_info_to_descriptors()
    molecules.append(mol1)

df = create_descriptor_df(molecules)
df.to_csv("rdkit_descriptors.csv", index=False)
