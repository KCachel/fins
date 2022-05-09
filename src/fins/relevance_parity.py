"""Metrics to assess the relevance parity fairness of a subset selection
    References
    ----------
    Kathleen Cachel and Elke Rundensteiner.
    "FINS Auditing Framework: Group Fairness for Subset Selections"
    in the proceedings of the AAAI/ACM conference on Artificial Intelligence,
    Ethics, and Society (AIES 2022)
"""


# Authors: Kathleen Cachel <kcachel@wpi.edu>
# License:


import numpy as np
import fins as fins


def relevance_parity(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups):
    """Compute the relevance parity.
    Parameters
    ----------
    pool_items : numpy array of shape = (n_items)
        The items in the pool (sorted by relevance score).
    pool_scores: numpy array of shape = (n_items)
        The scores of the items in the pool (sorted by relevance score).
    pool_groups: numpy array of shape = (n_items)
        The group identity of the items in the pool (corresponding to order of items in pool_items).
    subset_items : numpy array of shape = (n_items)
        The items in the subset(sorted by relevance score).
    subset_scores : numpy array of shape = (n_items)
        The scores of the items in the subset(sorted by relevance score).
    subset_groups: numpy array of shape = (n_items)
        The group identity of the items in the subset (corresponding to order of items in subset_items).
    Returns
    ----------
    RselectRt: numpy array of shape = (n_groups)
        The relevance select rate
    rp_val: float
        Statistical Parity fairness.
    Examples
    --------
    --------
    >>> pool_items = np.asarray([1,2,3,4])
    >>> pool_groups = np.asarray([0, 0, 1, 1])
    >>> pool_scores = np.asarray([100, 85, 54, 12])
    >>> subset_items = np.asarray([1,4])
    >>> subset_groups = np.asarray([0, 1])
    >>> subset_scores = np.asarray([100,12])
    >>> relevance_parity(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups)
    array([0.01081081, 0.03030303]), 0.3567567567567568
    """


    fins.check_pool_subset_groups(pool_items, pool_scores, pool_groups, subset_items, subset_scores,
                                 subset_groups)  # error handling
    unique_grps = np.unique(pool_groups)
    num_unique_grps = unique_grps.shape[0]
    RselectRt  = np.full((num_unique_grps,), -np.Inf)
    for grp in unique_grps:
        pool_mask = pool_groups == grp
        subset_mask = subset_groups == grp
        num_grp_items_in_subset = np.count_nonzero(subset_mask)
        num_grp_items_in_pool = np.count_nonzero(pool_mask)
        sum_scores_grp_in_pool = np.sum(pool_scores[pool_mask])
        avg_grp_score = sum_scores_grp_in_pool/ num_grp_items_in_pool
        RselectRt[grp] = num_grp_items_in_subset / avg_grp_score

    min_group_selection_prop  = np.min(RselectRt)
    max_group_selection_prop = np.max(RselectRt)
    rp_val = min_group_selection_prop / max_group_selection_prop
    return RselectRt, rp_val
