import pandas as pd
from scipy.stats import kstest
import matplotlib.pyplot as plt
import seaborn as sns

from typing import List, Optional


class FragmentAnalysis:
    def __init__(self, data: pd.DataFrame, fragment_columns: List, feature_columns: List):
        self.feature_columns = feature_columns
        self.fragment_columns = fragment_columns
        self.data = data

        self.feature_explanation = {} # key = feature

    def visualise_feature(self, feature: str, fragment: str):
        with_, without = self._get_with_and_without_fragment_datasets(feature, fragment)

        plt.figure(figsize=(15, 8))

        ax = sns.distplot(without, label="rest")
        ax = sns.distplot(with_, color="black", label="selected")
        return plt

    def _get_with_and_without_fragment_datasets(self, feature: str, fragment: str):
        with_ = self.data[self.data[fragment] != 0][feature]
        without = self.data[~(self.data[fragment] != 0)][feature]
        return with_, without

    def scan_fragments_for_explainability(self, feature: str):
        deltas = []
        for fragment in self.fragment_columns:
            with_, without = self._get_with_and_without_fragment_datasets(feature, fragment)
            if len(with_) == 0 or len(without)==0:
                continue
            compared = with_.describe().to_frame(f"with {fragment}")
            compared[f"without {fragment}"] = without.describe()
            compared[fragment] = abs(compared[f"with {fragment}"] - compared[f"without {fragment}"])

            stat = kstest(with_, without, alternative="two-sided")

            delta = compared[[fragment]].transpose()
            delta["K-S test p value"] = stat.pvalue
            delta["count"] = len(with_)
            mean_change = compared[f"with {fragment}"]["mean"] - compared[f"without {fragment}"]["mean"]
            median_change = compared[f"with {fragment}"]["50%"] - compared[f"without {fragment}"]["50%"]
            delta["mean change direction"] = "positive" if mean_change > 0 else "negative"
            delta["median change direction"] = "positive" if median_change > 0 else "negative"

            deltas.append(delta)

        explainability = pd.concat(deltas)
        explainability.rename(columns={"50%": "median"}, inplace=True)

        return explainability

    def find_fragment_explanation(self, feature: str, metric: Optional[str] = "median", n_frags: Optional[int]=10):

        explainability = self.scan_fragments_for_explainability(feature)
        explainability = explainability.sort_values(metric, ascending=False).reset_index()

        explanation = explainability.head(n_frags)

        return explanation[["index", metric, f"{metric} change direction"]]
