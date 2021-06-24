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
        """ Brace yourself this will take a while for now  In development..."""
        sim = SimilarityAnalysis(self.data, self.mol_column_name, self.target_column_name, features=[x, y])
        neighbors_list = []

        for i, mol in enumerate(self.data[self.mol_column_name]):
            neighbors = sim.find_similar_molecules(mol, n)
            print(f"Calculated {i}/{len(self.data[self.mol_column_name])}")
            neighbors_list.append(neighbors)

        neighbors_all = pd.concat(neighbors_list)
        neighbors_stats = neighbors_all[["molecule", "distance"]].groupby("molecule").mean()

        return neighbors_stats

    def plot_fragments_present(self, fragment_label: AnyStr, marker: AnyStr = None) -> object:
        """
        Plots fragment representation of the whole training dataset
        """
        # TODO REFACTOR WHOLE FUNCTION

        cols = [col for col in self.data.columns if fragment_label in col]
        all_fragments = self.data[cols]
        all_fragments = all_fragments.where(all_fragments[cols] == 0, 1)  # Converts to binary system.
        # all_fragments /= all_fragments # Converts to binary system.
        fragments = all_fragments[cols].sum().sort_values(ascending=False)
        fig, ax = plt.subplots(figsize=(20, 12))
        sns.set(style="darkgrid")

        if marker:
            mark = self.data[self.data[self.mol_column_name] == marker][fragments.index.to_list()]
            mark /= mark
            mark.replace(0, np.nan, inplace=True)
            line_widths = mark * 1.5
            bar1 = sns.barplot(x=fragments.index.to_list(), y=fragments.values, edgecolor='black',
                               linewidth=line_widths.values[0])
            plt.plot(mark.columns.values, mark.values[0], markersize=12, marker='x')
        else:
            bar1 = sns.barplot(x=fragments.index.to_list(), y=fragments.values)

        loc, labels = plt.xticks()
        bar1.set_xticklabels(labels, size=8, rotation=90)
        bar1.set_xlabel("Fragments")
        bar1.set_ylabel("Count")
        bar1.set_xlim(-1, len(cols))
        plt.tight_layout()
        return fig, ax

    def plot_feature_plane(self, feature_x: AnyStr, feature_y: AnyStr, marker: AnyStr = None) -> None:
        """
        Plots the feature plane for a given set of features.
        """
        # TODO REFACTOR WHOLE FUNCTION
        if self.data[feature_x].dtype not in ['float', 'float64', 'int', 'int64']:
            raise ValueError(f"{feature_x} is not a continous data feature please supply numerical features.")

        if self.data[feature_y].dtype not in ['float', 'float64', 'int', 'int64']:
            raise ValueError(f"{feature_y} is not a continous data feature please supply numerical features.")

        xmin, xmax = data[feature_x].min(), data[feature_x].max()
        ymin, ymax = data[feature_y].min(), data[feature_y].max()

        joint_plot = sns.jointplot(data=self.data, x=feature_x, y=feature_y,
                                   marginal_kws=dict(bins=100),
                                   xlim=(xmin, xmax), ylim=(ymin, ymax))

        joint_plot.ax_joint.cla()
        plt.sca(joint_plot.ax_joint)  # strip out the joint keeping the marginals

        plt.hist2d(data=self.data, x=feature_x, y=feature_y, bins=(100, 100), cmap="viridis",
                   range=[[xmin, xmax], [ymin, ymax]])
        if marker:
            mark = self.data[self.data[self.mol_column_name] == marker]
            plt.plot([mark[feature_x]], [mark[feature_y]], 'x', color='red', label=marker)
            plt.legend()
        plt.xlabel(feature_x)
        plt.ylabel(feature_y)


if __name__ == "__main__":
    data = pd.read_csv("../exploration/rdkit_descriptors.csv").dropna()
    features = [col for col in data.columns if 'fr' not in col]
    chem = ChemicalSpace(data=data, mol_column_name='SMILES', target_column_name='MP')
    a = chem.plot_feature_plane('qed', "MinAbsPartialCharge", marker="ClC1=CC(Cl)=C2N=CC=C(Br)C2=C1")
    b, c = chem.plot_fragments_present('fr', marker="ClC1=CC(Cl)=C2N=CC=C(Br)C2=C1")
