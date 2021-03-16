import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

from typing import Optional, Sequence
# conda install --file requirements.txt -c rdkit

class MoleculeAnalysis:
    def __init__(self, data: pd.DataFrame, features: Optional[Sequence[str]]=None):
        self.data = data
        self.feature_percentiles = pd.DataFrame()
        self.features = features

    def get_percentile(self, feature, mol):
        mol_data = self.data[self.data["SMILES"]==mol]
        mol_val = mol_data[feature].values[0]
        percentile = stats.percentileofscore(self.data[feature], mol_val)
        return percentile

    def get_percentiles_all_features(self, mol: str, feature_subset=None):
        if feature_subset:
            features = feature_subset
        else:
            features = self.features

        if mol not in self.feature_percentiles.columns:
            percentiles = []
            for feature in features:
                p = self.get_percentile(feature, mol)
                percentiles.append({"feature": feature, mol: p})
            percentiles = pd.DataFrame(percentiles)
            if self.feature_percentiles.empty:
                self.feature_percentiles = percentiles
            else:
                self.feature_percentiles = pd.merge(self.feature_percentiles, percentiles, on="feature")

        return self.feature_percentiles[["feature", mol]]

    def show_mol_in_dist(self, mol, feature, kde=False, hist=True):
        print(f'mol property value: {self.data[self.data["SMILES"] == mol][feature].values[0]}')
        plt.figure(figsize=(15, 8))
        ax = sns.distplot(self.data[feature], label=feature, kde=kde, hist=hist)
        plt.axvline(self.data[self.data["SMILES"] == mol][feature].values[0], color="black")
        return plt

    def show_mol_in_scatter(self, mol, x, y=None):
        if not y:
            y = "MP"

        plt.figure(figsize=(15, 8))
        ax = sns.scatterplot(x=self.data[x], y=self.data[y])
        plt.axvline(self.data[self.data["SMILES"] == mol][x].values[0], color="black")
        plt.axhline(self.data[self.data["SMILES"] == mol][y].values[0], color="black")

        return plt

    def select_features_from_percentiles(self, cut_off=0.15, tails="both", mols=None):
        if not mols:
            mols = [mol for mol in self.feature_percentiles if mol != "feature"]

    def draw_mol(self, mol):
        m = Chem.MolFromSmiles(mol)
        pic = Draw.MolToImage(m)
        return pic


def random():
    data = pd.read_csv("rdkit_descriptors.csv")
    data = data[data["qed"] != -666]

    frag_cols = []
    for col in list(data.columns):
        if col[:3] == "fr_":
            frag_cols.append(col)

    features = []
    for col in data.columns:
        if col not in frag_cols and col not in ["SMILES", "MP"]:
            features.append(col)

    MolAnalyser = MoleculeAnalysis(data, features)
    mol = "NC1=C2C=C(Cl)N=CC2=CC=C1"
    feature = "MolLogP"
    MolAnalyser.get_percentiles_all_features(mol)

    MolAnalyser.show_mol_in_dist(mol, feature)

    MolAnalyser.show_mol_in_scatter(mol, feature)

    MolAnalyser.draw_mol(mol)