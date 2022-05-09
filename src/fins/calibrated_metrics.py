"""Metrics to assess the calibrated balance fairness and the calibrated parity
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



def calibrated_parity(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups, lb_bin, ub_bin):
    """Compute the calibrated parity.
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
    lb_bin: numpy array of shape = (n_bins)
        The lower bound scores for each bin (bin is greater than or equal to lower bound).
    ub_bin: numpy array of shape = (n_bins)
        The upper bound scores for each bin (bin is less than upper bound).
    Returns
    ----------
    bin_group_selection_proportions: numpy array of shape = (n_bins,n_groups)
        The proportion of each group selected into the subset from the bin
    dp_val: float
        Calibrated parity value.
    Examples
    --------
    --------
    >>> pool_items = np.asarray([1,2,3,4])
    >>> pool_scores = np.asarray([100, 85, 54, 12])
    >>> pool_groups = np.asarray([0, 0, 1, 1])
    >>> subset_items = np.asarray([1,4])
    >>> subset_scores = np.asarray([100,12])
    >>> subset_groups = np.asarray([0, 1])
    >>> lb_bin = np.asarray([0, 50])
    >>> ub_bin = np.asarray([49, 100])

    >>> calibrated_parity(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups, lb_bin, ub_bin)
    (array([[0. , 1. ],
       [0.5, 0. ]]), array([0., 0.]), 0.0)
    """

    fins.check_pool_subset_groups(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups) #error handling
    n_bins = lb_bin.shape[0]
    unique_grps = np.unique(pool_groups)
    num_unique_grps = unique_grps.shape[0]
    bin_group_selectr = np.full((n_bins, num_unique_grps), -np.Inf)

    for bin_i in range(0,n_bins):
        lb = lb_bin[bin_i]
        ub = ub_bin[bin_i]
        greaterthanequal_lb_pool = pool_scores > lb
        lessthan_ub_pool = pool_scores <= ub
        bin_mask_pool = np.bitwise_and(greaterthanequal_lb_pool,lessthan_ub_pool)
        greaterthanequal_lb_subset = subset_scores > lb
        lessthan_ub_subset = subset_scores <= ub
        bin_mask_subset = np.bitwise_and(greaterthanequal_lb_subset, lessthan_ub_subset)
        bin_pool_items = pool_items[bin_mask_pool]
        bin_pool_groups = pool_groups[bin_mask_pool]
        bin_subset_items = subset_items[bin_mask_subset]
        bin_subset_groups = subset_groups[bin_mask_subset]
        for grp in unique_grps:
            grp_bin_pool_mask = bin_pool_groups == grp
            grp_bin_pool_items = bin_pool_items[grp_bin_pool_mask]
            num_grp_bin_pool_items = np.count_nonzero(grp_bin_pool_items)
            grp_bin_subset_mask = bin_subset_groups == grp
            grp_bin_subset_items = bin_subset_items[grp_bin_subset_mask]
            num_grp_bin_subset_items = np.count_nonzero(grp_bin_subset_items)
            if num_grp_bin_pool_items == 0:
                bin_group_selectr[bin_i, grp] = 0.0
            else:
                bin_group_selectr[bin_i, grp] = num_grp_bin_subset_items / num_grp_bin_pool_items

    max_prop_each_bin = np.max(bin_group_selectr, axis = 1)
    min_prop_each_bin = np.min(bin_group_selectr, axis = 1)


    if np.all(min_prop_each_bin == max_prop_each_bin):
        cp_val = 1 #totally fair since max = min in all bins
    else:
        different_bin_selection_rates_mask = min_prop_each_bin != max_prop_each_bin
        cp_val = np.min(min_prop_each_bin[different_bin_selection_rates_mask] / max_prop_each_bin[different_bin_selection_rates_mask])
    return bin_group_selectr, cp_val


def calibrated_balance(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups, lb_bin, ub_bin):
    """Compute the calibrated balance.
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
    lb_bin: numpy array of shape = (n_bins)
        The lower bound scores for each bin (bin is greater than or equal to lower bound).
    ub_bin: numpy array of shape = (n_bins)
        The upper bound scores for each bin (bin is less than upper bound).
    Returns
    ----------
    bin_group_proportions: numpy array of shape = (n_bins,n_groups)
        The proportion of each group selected into the subset from the bin
    db_val: float
        Distributed parity value.
    Examples
    --------
    --------
    >>> pool_items = np.asarray([1,2,3,4])
    >>> pool_scores = np.asarray([100, 85, 54, 12])
    >>> pool_groups = np.asarray([0, 0, 1, 1])
    >>> subset_items = np.asarray([2,4])
    >>> subset_scores = np.asarray([85,12])
    >>> subset_groups = np.asarray([0, 1])
    >>> lb_bin = np.asarray([0, 87])
    >>> ub_bin = np.asarray([86, 100])

    >>> calibrated_balance(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups, lb_bin, ub_bin)
    (array([[0.5, 0.5],
       [0. , 0. ]]), 1)
    """
    fins.check_pool_subset_groups(pool_items, pool_scores, pool_groups, subset_items, subset_scores,
                                 subset_groups)  # error handling
    n_bins = lb_bin.shape[0]
    unique_grps = np.unique(pool_groups)
    num_unique_grps = unique_grps.shape[0]
    bin_group_proportions = np.full((n_bins, num_unique_grps), -np.Inf)

    for bin_i in range(0,n_bins):
        lb = lb_bin[bin_i]
        ub = ub_bin[bin_i]
        greaterthanequal_lb_pool = pool_scores > lb
        lessthan_ub_pool = pool_scores <= ub
        bin_mask_pool = np.bitwise_and(greaterthanequal_lb_pool,lessthan_ub_pool)
        greaterthanequal_lb_subset = subset_scores > lb
        lessthan_ub_subset = subset_scores <= ub
        bin_mask_subset = np.bitwise_and(greaterthanequal_lb_subset, lessthan_ub_subset)
        bin_pool_items = pool_items[bin_mask_pool]
        bin_pool_groups = pool_groups[bin_mask_pool]
        bin_subset_items = subset_items[bin_mask_subset]
        bin_subset_groups = subset_groups[bin_mask_subset]
        for grp in unique_grps:
            grp_bin_pool_mask = bin_pool_groups == grp
            grp_bin_pool_items = bin_pool_items[grp_bin_pool_mask]
            num_grp_bin_pool_items = np.count_nonzero(grp_bin_pool_items)
            grp_bin_subset_mask = bin_subset_groups == grp
            grp_bin_subset_items = bin_subset_items[grp_bin_subset_mask]
            num_grp_bin_subset_items = np.count_nonzero(grp_bin_subset_items)
            num_subset_items = np.count_nonzero(subset_items)
            if num_grp_bin_pool_items == 0:
                bin_group_proportions[bin_i, grp] = 0.0
            else:
                bin_group_proportions[bin_i, grp] = num_grp_bin_subset_items / num_subset_items

    max_props_each_bin = np.max(bin_group_proportions, axis = 1)
    min_props_each_bin = np.min(bin_group_proportions, axis = 1)




    if np.all(min_props_each_bin == max_props_each_bin):
        cb_val = 1 #totally fair since max = min in all bins
    else:
        different_bin_selection_rates_mask = min_props_each_bin != max_props_each_bin
        cb_val = np.min(min_props_each_bin[different_bin_selection_rates_mask] / max_props_each_bin[different_bin_selection_rates_mask])
    return bin_group_proportions, cb_val

