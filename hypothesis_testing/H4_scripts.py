from sklearn.neighbors import NearestNeighbors
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw


class SimilarityAnalysis:
    def __init__(self, data: pd.DataFrame, mol_column_name, target_column_name, features=None):
        self.data = data
        if features is None:
            features = [col for col in data.columns if col not in [mol_column_name, target_column_name]]
        self.mol_column = mol_column_name
        self.target_column = target_column_name
        self.key_columns = [mol_column_name, target_column_name]
        self.features = features

        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(data[features])
        self.nbrs = nbrs

    def _process_neighbors(self, distances, indices, smiles, return_self=False):
        d = pd.DataFrame(distances.T, columns=["distance"])
        d["id"] = indices.T
        d = d.merge(self.data[self.key_columns], left_on="id", right_index=True)
        d.rename(columns={self.mol_column: "similar"}, inplace=True)
        d["molecule"] = smiles
        if not return_self:
            d = d[d["distance"] != 0]
        return d

    def find_similar_molecules(self, smiles, n_similar=5):
        mol_data = self.data[self.data[self.mol_column] == smiles][self.features]
        distances, indices = self.nbrs.kneighbors(mol_data, n_neighbors=n_similar, return_distance=True)
        neighbors = self._process_neighbors(distances, indices, smiles)
        return neighbors

    def find_molecules_in_radius(self, smiles, radius=1):
        # should probably make this accept multiple mols at once
        mol_data = self.data[self.data[self.mol_column] == smiles][self.features]
        distances, indices = self.nbrs.radius_neighbors(mol_data, radius=radius, return_distance=True)
        distances, indices = distances[0], indices[0] # weird array in a array situation
        neighbors = self._process_neighbors(distances, indices, smiles)

        return neighbors

    @staticmethod
    def generate_visualisation(df):
        mol_dict = df.groupby('molecule')["similar"].apply(list).to_dict()

        all_visuals = []
        for mol, similars in mol_dict.items():
            m = Chem.MolFromSmiles(mol)
            pic = Draw.MolToImage(m)
            ms = []
            for similar in similars:
                ms.append(Chem.MolFromSmiles(similar))

            pics = Draw.MolsToGridImage(ms)
            all_visuals.append((pic, pics))

        return all_visuals

    def analyse_neighbors(self, smiles, radius=1):
        if isinstance(smiles, str):
            smiles = [smiles]

        neighbors_list = []
        for mol in smiles:
            neighbors = self.find_molecules_in_radius(mol, radius)
            neighbors_list.append(neighbors)

        neighbors_all = pd.concat(neighbors_list)

        neighbors_stats = neighbors_all[["molecule", self.target_column]].groupby("molecule").describe()

        return neighbors_stats

