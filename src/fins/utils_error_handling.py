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

# Function for error handling with group information
def check_pool_subset_groups(pool_items, pool_scores, pool_groups, subset_items, subset_scores, subset_groups):
    """
        Check whether input pool, subset, and protected attribute information is valid.

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
    Raise errors if found.
    """

    # error handling for types
    if not isinstance(pool_items, np.ndarray):
        raise TypeError("Input pool_items must be ndarray")

    if not isinstance(pool_scores, np.ndarray):
        raise TypeError("Input pool_scores must be ndarray")

    if not isinstance(pool_groups, np.ndarray):
        raise TypeError("Input pool_groups must be ndarray")

    if not isinstance(subset_items, np.ndarray):
        raise TypeError("Input subset_items must be ndarray")

    if not isinstance(subset_scores, np.ndarray):
        raise TypeError("Input subset_scores must be ndarray")

    if not isinstance(subset_groups, np.ndarray):
        raise TypeError("Input subset_groups must be ndarray")


    # error handling for inputs sizes and values
    if len(pool_items) <= 0:  # check size of input pool
        raise ValueError("Please input a valid pool")
    if len(subset_items) <= 0:  # check size of subset
        raise ValueError("Please input a valid subset")

    if len(pool_items) != len(pool_scores):  # check size of input pool and scores match
        raise ValueError("Please input a pool and scores of same size")

    if len(pool_items) != len(pool_groups):  # check size of input pool and groups match
        raise ValueError("Please input a pool and groups of same size")

    if len(subset_items) != len(subset_scores):  # check size of input subset and scores match
        raise ValueError("Please input a subset and scores of same size")

    if len(subset_items) != len(subset_groups):  # check size of input subset and groups match
        raise ValueError("Please input a subset and groups of same size")

    if len(set(pool_items)) != len(pool_items):  # check for repetition in input pool
        raise ValueError("Please input a pool with unique items")

    if len(set(subset_items)) != len(subset_items):  # check for repetition in subset
        raise ValueError("Please input a pool with unique items")

    if (np.min(np.unique(pool_groups)) != 0) or (np.max(np.unique(pool_groups)) != len(
            np.unique(pool_groups)) - 1):  # check for consecutive group encoding in pool
        raise ValueError("Please represent pool groups via consecutive integers 0,1,2, etc...")


# Function for error handling with group information
def check_pool_subset_items_groups(pool_items, pool_groups, subset_items, subset_groups):
    """
        Check whether input pool, subset, and protected attribute information is valid.

        Parameters
    ----------
    pool_items : numpy array of shape = (n_items)
        The items in the pool (sorted by relevance score).
    pool_groups: numpy array of shape = (n_items)
        The group identity of the items in the pool (corresponding to order of items in pool_items).
    subset_items : numpy array of shape = (n_items)
        The items in the subset(sorted by relevance score).
    subset_groups: numpy array of shape = (n_items)
        The group identity of the items in the subset (corresponding to order of items in subset_items).

    Returns
    ----------
    Raise errors if found.
    """

    # error handling for types
    if not isinstance(pool_items, np.ndarray):
        raise TypeError("Input pool_items must be ndarray")

    if not isinstance(pool_groups, np.ndarray):
        raise TypeError("Input pool_groups must be ndarray")

    if not isinstance(subset_items, np.ndarray):
        raise TypeError("Input subset_items must be ndarray")


    if not isinstance(subset_groups, np.ndarray):
        raise TypeError("Input subset_groups must be ndarray")


    # error handling for inputs sizes and values
    if len(pool_items) <= 0:  # check size of input pool
        raise ValueError("Please input a valid pool")
    if len(subset_items) <= 0:  # check size of subset
        raise ValueError("Please input a valid subset")

    if len(pool_items) != len(pool_groups):  # check size of input pool and groups match
        raise ValueError("Please input a pool and groups of same size")

    if len(subset_items) != len(subset_groups):  # check size of input subset and groups match
        raise ValueError("Please input a subset and groups of same size")

    if len(set(pool_items)) != len(pool_items):  # check for repetition in input pool
        raise ValueError("Please input a pool with unique items")

    if len(set(subset_items)) != len(subset_items):  # check for repetition in subset
        raise ValueError("Please input a pool with unique items")

    if (np.min(np.unique(pool_groups)) != 0) or (np.max(np.unique(pool_groups)) != len(np.unique(pool_groups))- 1): #check for consecutive group encoding in pool
        raise  ValueError("Please represent pool groups via consecutive integers")


# Function for error handling with group information
def check_subset_items_groups(pool_groups, subset_items, subset_groups):
    """
        Check whether input pool gorups, subset, and protected attribute information is valid.

        Parameters
    ----------
    pool_groups: numpy array of shape = (n_items)
        The group identity of the items in the pool (corresponding to order of items in pool_items).
    subset_items : numpy array of shape = (n_items)
        The items in the subset(sorted by relevance score).
    subset_groups: numpy array of shape = (n_items)
        The group identity of the items in the subset (corresponding to order of items in subset_items).

    Returns
    ----------
    Raise errors if found.
    """

    # error handling for types
    if not isinstance(pool_groups, np.ndarray):
        raise TypeError("Input pool_groups must be ndarray")

    if not isinstance(subset_items, np.ndarray):
        raise TypeError("Input subset_items must be ndarray")

    if not isinstance(subset_groups, np.ndarray):
        raise TypeError("Input subset_groups must be ndarray")

    # error handling for inputs sizes and values
    if len(subset_items) <= 0:  # check size of subset
        raise ValueError("Please input a valid subset")

    if len(subset_items) != len(subset_groups):  # check size of input subset and groups match
        raise ValueError("Please input a subset and groups of same size")

    if len(set(subset_items)) != len(subset_items):  # check for repetition in subset
        raise ValueError("Please input a pool with unique items")

    if (np.min(np.unique(pool_groups)) != 0) or (np.max(np.unique(pool_groups)) != len(
            np.unique(pool_groups)) - 1):  # check for consecutive group encoding in pool
        raise ValueError("Please represent pool groups via consecutive integers")