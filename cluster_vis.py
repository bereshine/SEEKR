"""
Visualization for clusters

@author: Chris Horning

"""

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import optimal_leaf_ordering
from scipy.cluster.hierarchy import leaves_list

def cluster_kmers(kmers):
    """

    :param kmers: numpy.ndarray of kmer counts
    :return: ndarray linkage matrix
    """

    #see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html

    #Note: We are splitting up the operations of linkage for clarity and benchmarking, these could be combined into
    #one call to linkage

    #compute the pairwise distance between sequences.
    y = pdist(kmers, metric='correlation')
    Z = linkage(y, method='ward', optimal_ordering=False)
    Z_ordered = optimal_leaf_ordering(Z, y)

    return Z_ordered

def get_ordering(Z):

    leaves = leaves_list(Z)
    return leaves

