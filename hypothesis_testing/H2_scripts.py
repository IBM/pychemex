import pandas as pd

from typing import List, Optional


class FragmentAnalysis:
    def __init__(self, data: pd.DataFrame, fragment_columns: List, feature_columns: List):
        self.feature_columns = feature_columns
        self.fragment_columns = fragment_columns
        self.data = data

    def visualise_feature(self, feature: str, fragment: str):
        pass

    def scan_features(self, fragment: str):
        pass

    def find_explainable_features(self, fragment: str):
        pass

"""
with_ = data[data[fragment]!=0][variable]
without = data[~(data[fragment]!=0)][variable]

compared = with_.describe().to_frame(f"with {fragment}")
compared[f"without {fragment}"] = without.describe()
compared[fragment] = abs(compared[f"with {fragment}"] - compared[f"without {fragment}"])

plt.figure(figsize=(15, 8))

ax = sns.distplot(without, label="rest")
ax = sns.distplot(with_, color="black", label="selected")

"1. go thru features and get the difference in mean for have features vs ~(have feature)\n",
    "    1. delta mean\n",
    "    2. delta median\n",
    "    3. Mann Whitney U test - give warning for unbalanced datasets... maybe not use at all since it's likely to be unbalanced\n",
    "    4. prob that subset comes from population\n",
    "2. select a feature\n",
    "3. plot features vs ~(have feature)\n",
    "4. repeat by scanning thru fragments?\n",
    "5. select X features that are significant for fragment Y\n",
    "6. think of how to handle categorical features (binary, ordered, unordered)"
"""