# Forms part of H10 - What is your chemical space
# https://github.com/IBM/pychemex/issues/19

from typing import AnyStr

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from hypothesis_testing.H4_scripts import SimilarityAnalysis


class ChemicalSpace:
    """
    Plots the chemical space of feature pairs.
    """

    def __init__(self, data: pd.DataFrame, mol_column_name: AnyStr, target_column_name: AnyStr):
        """

        :param data: Dataframe containing all features and target
        :param mol_column_name: Name of molecule identifier column
        :param target_column_name: Name of target column
        """
        self.data = data
        self.mol_column_name = mol_column_name
        self.target_column_name = target_column_name

    def _compute_neighbors(self, x, y, n):
        """ Brace yourself this will take a while for now """
        sim = SimilarityAnalysis(self.data, self.mol_column_name, self.target_column_name, jobs=-1, features=[x, y])
        neighbors_list = []

        for mol in self.data[self.mol_column_name]:
            neighbors = sim.find_similar_molecules(mol, n)
            neighbors_list.append(neighbors)

        neighbors_all = pd.concat(neighbors_list)
        neighbors_stats = neighbors_all[["molecule", "distance"]].groupby("molecule").mean()

        return neighbors_stats

    def plot_feature_plane(self, feature_x: AnyStr, feature_y: AnyStr):
        """
        Plots the feature plane
        """

        if self.data[feature_x].dtype not in ['float', 'float64', 'int', 'int64']:
            raise ValueError(f"{feature_x} is not a continous data feature please supply numerical features.")

        if self.data[feature_y].dtype not in ['float', 'float64', 'int', 'int64']:
            raise ValueError(f"{feature_y} is not a continous data feature please supply numerical features.")

        xmin, xmax = data[feature_x].min(), data[feature_x].max()
        ymin, ymax = data[feature_y].min(), data[feature_y].max()


        all_values = self._compute_neighbors(feature_x, feature_y,5)
        self.data = self.data.merge(all_values, left_on="SMILES", right_index=True)
        self.data["Density_coeff"] = 1/ (self.data['distance']+ 1e-18)
        self.data["Density_coeff"] = self.data["Density_coeff"].fillna(0)
        ax = sns.jointplot(data=self.data, x=feature_x, y=feature_y,
                           marginal_kws=dict(bins=100),
                           xlim=(0, xmax), ylim=(0, ymax))

        ax.ax_joint.cla()
        plt.sca(ax.ax_joint)  # strip out the joint keeping the marginals


        plt.tricontourf(self.data[feature_x], self.data[feature_y], self.data["Density_coeff"])
        # plt.hist2d(data=self.data, x=feature_x, y=feature_y, bins=(100, 100), cmap=cm.viridis,
        #            range=[[0, xmax], [0, ymax]])
        plt.colorbar()
        plt.xlim()
        plt.xlabel(feature_x)
        plt.ylabel(feature_y)
        return all_values


if __name__ == "__main__":
    data = pd.read_csv("../exploration/rdkit_descriptors.csv")
    features = [col for col in data.columns if 'fr' not in col]
    chem = ChemicalSpace(data=data[:10000], mol_column_name='SMILES', target_column_name='MP')
    a = chem.plot_feature_plane('TPSA', "MolWt")
