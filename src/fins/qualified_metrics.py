"""Metrics to assess the qualified balance fairness and the qualified parity
 fairness of a subset selection
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

def qualififed_parity(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups, q):
    """Compute qualified parity.
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
    q: float
        The relevance score for which items in the pool that have score >= q are "relevant".
    Returns
    ----------
    QselectRt: numpy array of shape = (n_groups)
        The proportion of each group selected into the subset from the qualified pool
    qp_val: float
        fairness.
    Examples
    --------
    --------
    >>> pool_items = np.asarray([1,2,3,4])
    >>> pool_scores = np.asarray([100, 85, 54, 12])
    >>> pool_groups = np.asarray([0, 0, 1, 1])
    >>> subset_items = np.asarray([2,4])
    >>> subset_scores = np.asarray([85,12])
    >>> subset_groups = np.asarray([0, 1])
    >>> q = 50

array([0.5, 0. ]), 0.0
    """

    fins.check_pool_subset_groups(pool_items, pool_scores, pool_groups, subset_items, subset_scores,
                                 subset_groups)  # error handling
    unique_grps = np.unique(pool_groups)
    num_unique_grps = unique_grps.shape[0]
    QselectRt  = np.full((num_unique_grps,), -np.Inf)
    qualified_mask = pool_scores >= q
    qualified_pool_groups =  pool_groups[qualified_mask]
    subset_qualified_mask = subset_scores >= q
    qualified_subset_groups = subset_groups[subset_qualified_mask]
    for grp in unique_grps:
        pool_mask = qualified_pool_groups == grp
        subset_mask = qualified_subset_groups == grp
        num_grp_items_in_subset = np.count_nonzero(subset_mask)
        num_grp_items_in_pool = np.count_nonzero(pool_mask)
        if num_grp_items_in_pool == 0:
            QselectRt[grp] = 0.0
        else:
            QselectRt[grp] = num_grp_items_in_subset /num_grp_items_in_pool

    min_group_selection_prop  = np.min(QselectRt)
    max_group_selection_prop = np.max(QselectRt)
    qp_val =  min_group_selection_prop / max_group_selection_prop



    return QselectRt, qp_val


def qualified_balance(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups, q):
    """Compute qualified balance.
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
    q: float
        The relevance score for which items in the pool that have score >= q are "relevant".
    Returns
    ----------
    QpropOfS: numpy array of shape = (n_groups)
        The proportion of each group selected into the subset from the qualified pool
    qb_val: float
        qualified balance fairness.
    Examples
    --------
    --------
    >>> pool_items = np.asarray([1,2,3,4])
    >>> pool_scores = np.asarray([100, 85, 54, 12])
    >>> pool_groups = np.asarray([0, 0, 1, 1])
    >>> subset_items = np.asarray([2,4])
    >>> subset_scores = np.asarray([85,12])
    >>> subset_groups = np.asarray([0, 1])
    >>> q = 50

array([0.5, 0. ]), 0.0
    """
    fins.check_pool_subset_groups(pool_items, pool_scores, pool_groups, subset_items, subset_scores,
                                 subset_groups)  # error handling
    unique_grps = np.unique(pool_groups)
    num_unique_grps = unique_grps.shape[0]
    QpropOfS  = np.full((num_unique_grps,), -np.Inf)
    qualified_mask = pool_scores >= q
    qualified_pool_groups =  pool_groups[qualified_mask]
    subset_qualified_mask = subset_scores >= q
    qualified_subset_groups = subset_groups[subset_qualified_mask]
    total_items_subset = subset_items.shape[0]
    for grp in unique_grps:
        pool_mask = qualified_pool_groups == grp
        subset_mask = qualified_subset_groups == grp
        num_grp_items_in_subset = np.count_nonzero(subset_mask)
        num_grp_items_in_pool = np.count_nonzero(pool_mask)
        if num_grp_items_in_pool == 0:
            QpropOfS[grp] = 0.0
        else:
            QpropOfS[grp] = num_grp_items_in_subset /total_items_subset

    min_group_selection_prop  = np.min(QpropOfS)
    max_group_selection_prop = np.max(QpropOfS)
    qb_val =  min_group_selection_prop / max_group_selection_prop



    return QpropOfS, qb_val

