# fins
Group fairness auditing methods for set selections



**fins** is a Python library that provides group fairness auditing metrics for a variety of subset selection problems. The package includes a suite of metrics.



## Table of Contents
1. [Basic installation instructions](#basic-installation-instructions)
2. [Quick start examples](#quick-start-examples)
3. [Metrics in fins](#mettics-in-find)
4. [Citing fins](#citing-fins)



## Basic installation instructions
1. (Optionally) create a virtual environment
```
python3 -m venv fins
source fins/bin/activate
```
2. Install via pip
```
pip install finsfairauditing
```
You can also install fins directly from source.
```
git clone https://github.com/Kcachel/fins.git
cd fins
pip install -r requirements.txt
```

## Quick start examples
Fins contains a suite of group fairness metrics for subset selection tasks.
For a sample auditing guide see [fins github repo](https://github.com/KCachel/fins).

To use fins in your code run

```py
from finsfairauditing import fins
```

To audit a subset. Create the larger pool and subset:

```py
import numpy as np
pool_items = np.asarray([1,2,3,4])
pool_groups = np.asarray([0, 0, 1, 1])
pool_scores = np.asarray([100, 85, 54, 12])
subset_items = np.asarray([1,4])
subset_groups = np.asarray([0, 1])
subset_scores = np.asarray([100,12])
```
Choose the relevant fairness metric. Then return the per group metric and the high-level metric.
```py
rselectrt, rp_val = fins.relevance_parity(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups)
```
## Metrics in fins
fins provides the following pre-defined group fairness metrics for subset selectionss:
- **Parity**: statistical parity (proportional presence) group fairness of the selected set. To audit if the selected set contains a proportional number of items from each group.
- **Balance**: equal presence group fairness of the selected set. To audit if the selected set contains an equal number of items from each group.
- **Conditioned Parity**: statistical parity (proportional presence) group fairness of the selected set for all items sharing an additional attribute value(e.g., statistical parity for interns vs. whole population). To audit if the selected set contains a proportional number of items from each group conditional on items sharing some additional value.
- **Conditioned Balance**: equal presence group fairness of the selected set for all items sharing an additional attribute value (e.g., balance for interns vs. whole population). To audit if the selected set contains an equal number of items from each group conditional on items sharing some additional value.
- **Qualified Parity**: statistical parity (proportional presence) group fairness of the selected set for items deemed qualified (i.e., score greater than or equal to  q). To audit if the selected set contain a proportional presence of qualified items from each group.
- **Qualified Balance**: equal presence group fairness of the selected set (i.e., score greater than or equal to  q). To audit if the selected set contain an equal number of qualified items from each group.
- **Calibrated Parity**: statistical parity (proportional presence) group fairness of the selected set from specified score bins. To audit if items with similiar scores are if items with similar scores are treated similarly (via proportional presence) regardless of group membership.
- **Calibrated Balance**: equal presence group fairness of the selected set.  To audit if items with similiar scores are if items with similar scores are treated similarly (via equal presence) regardless of group membership.
- **Relevancce Parity**: To audit if groups are represented proportional to their average score (i.e., score-based relevance).
- **Score Parity**: To audit if the group-total score of the selected set is proportional to the number of items per group in the set.
- **Score Balance**: To audit if each groups receive and equal share of the selected set's total score.

## Citing fins
If you uses the fins auditing package we encourage
you to cite our paper:
```
@inproceedings{cachel2022fins,
  title={FINS Auditing Framework: Group Fairness for Subset Selections},
  author={Cachel, Kathleen and Rundensteiner, Elke},
  booktitle={Proceedings of the 2021 AAAI/ACM Conference on AI, Ethics, and Society},
  year={2022}
}
```
