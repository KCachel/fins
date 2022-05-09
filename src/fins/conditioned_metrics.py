"""Metrics to assess the conditioned balance fairness and the calibrated parity
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


def cond_parity(pool_items, pool_groups, pool_attribute, subset_items, subset_groups, subset_attribute, cond_val):
    """Compute conditioned parity.
    Parameters
    ----------
    pool_items : numpy array of shape = (n_items)
        The items in the pool (sorted by relevance score).
    pool_groups: numpy array of shape = (n_items)
        The group identity of the items in the pool (corresponding to order of items in pool_items).
    pool_attribute: numpy array of shape = (n_items)
        The additional attribute of the items in the pool (corresponding to order of items in pool_items).
    subset_items : numpy array of shape = (n_items)
        The items in the subset(sorted by relevance score).
    subset_groups: numpy array of shape = (n_items)
        The group identity of the items in the subset (corresponding to order of items in subset_items).
    subset_attribute: numpy array of shape = (n_items)
        The additional attribute of the items in the subset (corresponding to order of items in subset_items).
    cond_val: value to condition on
    Returns
    ----------
    CselectRt: numpy array of shape = (n_groups)
        The proportion of each conditioned group selected into the subset
    csp_val: float
        Conditioned Statistical Parity fairness.
    Examples
    --------
    --------
    >>> pool_items = np.asarray([1,2,3,4])
    >>> pool_groups = np.asarray([0, 0, 1, 1])
    >>> pool_attribute = np.asarray([0, 1, 0, 1])
    >>> subset_items = np.asarray([1,4])
    >>> subset_groups = np.asarray([0, 1])
    >>> subset_attribute = np.asarray([0, 1])
    >>> csp_val = 1
    >>> cond_parity(pool_items, pool_groups, pool_attribute, subset_items, subset_groups, subset_attribute, cond_val)
    array([0., 1.]), 0.0
    """

    fins.check_pool_subset_items_groups(pool_items, pool_groups, subset_items, subset_groups) #error handling
    unique_grps = np.unique(pool_groups)
    num_unique_grps = unique_grps.shape[0]
    CselectRt  = np.full((num_unique_grps,), -np.Inf)
    #narrow the pool
    cond_pool_mask = pool_attribute == cond_val
    pool_groups = pool_groups[cond_pool_mask]
    #narrow the subset
    cond_subset_mask = subset_attribute == cond_val
    subset_groups = subset_groups[cond_subset_mask]
    for grp in unique_grps:
        pool_mask = pool_groups == grp
        subset_mask = subset_groups == grp
        num_grp_items_in_subset = np.count_nonzero(subset_mask)
        num_grp_items_in_pool = np.count_nonzero(pool_mask)
        CselectRt[grp] = num_grp_items_in_subset /num_grp_items_in_pool

    min_group_selection_prop  = np.min(CselectRt)
    max_group_selection_prop = np.max(CselectRt)
    csp_val = min_group_selection_prop / max_group_selection_prop
    return CselectRt, csp_val

def cond_balance(pool_groups, subset_items, subset_groups, subset_attribute, cond_val):
    """Compute the conditioned balance.
    Parameters
    ----------
    pool_groups: numpy array of shape = (n_items)
        The group identity of the items in the pool (corresponding to order of items in pool_items).
    subset_items : numpy array of shape = (n_items)
        The items in the subset(sorted by relevance score).
    subset_groups: numpy array of shape = (n_items)
        The group identity of the items in the subset (corresponding to order of items in subset_items).
    subset_attribute: numpy array of shape = (n_items)
        The additional attribute of the items in the subset (corresponding to order of items in subset_items).
    cond_val: value to condition on
    Returns
    ----------
    CpropOfS: numpy array of shape = (n_groups)
        The proportion of the subset each group is
    cbal_val: float
        The balance error
    Examples
    --------
    --------
    >>> pool_groups = np.asarray([0, 0, 1, 1])
    >>> subset_items = np.asarray([1,4])
    >>> subset_groups = np.asarray([0, 1])
    >>> subset_attribute = np.asarray([1, 1])
    >>> cond_val = 1
    >>> cond_balance(pool_groups, subset_items, subset_groups, subset_attribute, cond_val)
    [0.5 0.5] 0.0
    """

    fins.check_subset_items_groups(pool_groups, subset_items, subset_groups)  # error handling
    unique_grps = np.unique(pool_groups)
    num_unique_grps = unique_grps.shape[0]
    CpropOfS  = np.full((num_unique_grps,), -np.Inf)

    # narrow the subset
    cond_subset_mask = subset_attribute == cond_val
    subset_items = subset_items[cond_subset_mask]
    subset_groups = subset_groups[cond_subset_mask]
    total_items_subset = subset_items.shape[0]
    for grp in unique_grps:
        subset_mask = subset_groups == grp
        num_grp_items_in_subset = np.count_nonzero(subset_mask)
        CpropOfS[grp] = num_grp_items_in_subset /total_items_subset

    min_group_proportion_subset  = np.min(CpropOfS)
    max_group_proportion_subset = np.max(CpropOfS)
    cbal_val =  min_group_proportion_subset / max_group_proportion_subset



    return CpropOfS, cbal_val



