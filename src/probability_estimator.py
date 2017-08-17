"""
This file contains a probability estimator method. Given a phenotype-genotype
pair within the context of a graph, and some statistics for different types
of path linking connected and disconnected pairs within the graph, it estimates
the probability that the pair is directly connected by a gene.
This file also contains a method for extracting the aforementioned path type
statistics from some data files, representing the alternative paths that link
together some phenotype-genotype pairs from a subset of the graph.
"""

from __future__ import division
from multiprocessing import Pool, cpu_count
from scipy.stats import norm
from graph_methods import get_phenotype_genotype_path_code, find_all_paths
import numpy as np
import operator
import time


def get_stats_from_alternative_paths(total_results_path, type_index_path, list_elems_path):
    """
    Loads from some files a set of alternative paths between phenotype-genotype
    pairs from a graph, return some statistics.
    Intended for collecting statistics from pairs that are known to be
    either all connected or all disconnected.

    Args:
        -total_results_path: Location of the file where paths are stored
            +Type: str
            +Pickle type: 2D numpy.array. Rows are phenotype-genotype pairs.
                Columns are path types
        - list_elems: Location of the file where list of computed pairs are stored.
                    Serves as row index.
            +Type: str
            +Pickle type: list[(phenotype_id, genotype_id)]
        -type_index: Location of the file where the dictionary of path type
                    and column index is stored.
            +Type: str
            +Pickle type: dict{path_type, index}

    Returns:
        -statistics: Mean and standard deviation for each path type
            +Type: list[(type, mean, stdDev)]
    """
    # Load the three structures
    # Type: 2D numpy.array. Rows are phenotype-genotype pairs.
    total_results = pickle.load(open(total_results_path, 'rb'))
    # Type: dict{path_type, index}
    type_index = pickle.load(open(type_index_path, 'rb'))
    # Type: list[(phenotype_id ,genotype_id)]
    list_elems = pickle.load(open(list_elems_path, 'rb'))

    # Compute mean and stddev for each path type
    means = np.mean(total_results, axis=0)
    stddevs = np.std(total_results, axis=0)
    # Sort path types by index (i.e., dictonary value)
    sorted_types = sorted(type_index.items(), key=operator.itemgetter(1))
    # Zip the statistics as a single function
    statistics = zip([x[0] for x in sorted_types], means, stddevs)
    return statistics


def estimate_probabilities(graph, common_genes, pair, statistics_conn, statistics_disc):
    """
    For a given pair of phenotype source and genotype target in a graph,
    use statistics for connected and disconnected pairs to estimate the
    probability that the pair is connected/disconnected.
    The probability estimator is based on normal distributions, with their
    means and standard deviations set to match those obtained for each type of
    path in the connected and disconnected cases.

    Args:
        -G: Graph where the paths are being computed.
            +Type: igraph graph
        -common_genes: List of genes linked with source phenotype that are ignored by G
            +Type: list[str]
        -pair: The source phenotype and target genotype ids.
            +Type: tuple(phenotype_id, genotype_id)
        -statistics_conn: Mean and standard deviation for each path type (connected)
            +Type: list[(type, mean, stdDev)]
        -statistics_disc: Mean and standard deviation for each path type (disconnected)
            +Type: list[(type, mean, stdDev)]

    Returns:
        -probability_conn: The probability that the pair is connected
            +Type: float
        -probability_disc: The probability that the pair is disconnected
            +Type: float
    """
    # Find all paths of length 4 between the source and target
    paths = find_all_paths(graph, graph.vs.find(pair[0]).index, graph.vs.find(pair[1]).index, maxlen=4)
    # Compute the type and frequency of all found paths
    found_paths = {}
    for current_path in paths:
        current_code = get_phenotype_genotype_path_code(graph, current_path, pair[0], common_genes, pair[1])
        if current_code in found_paths:
            found_paths[current_code] += 1
        else:
            found_paths[current_code] = 1
    # Create a set of all path codes, found and/or in the statistics.
    found_path_codes = set(found_paths.keys())
    all_path_codes = found_path_codes | set(x[0] for x in statistics_conn) | set(x[0] for x in statistics_disc)
    # Compute the probabilities for each path type
    probability_conn = 1  # Product of all connected probabilities
    probability_disc = 1  # Product of all disconnected probabilities
    path_num = 0
    for code in all_path_codes:
        # If paths of that path type were found for the pair, get how many
        if code in found_path_codes:
            path_num = found_paths[code]
        else:
            path_num = 0
        # Check if the type is found in the connected statistics.
        # If not found, we assume mean = 0, and stdDev = 1.
        mean = 0
        std = 1
        #If path statistics in connected data, load it
        for stat in statistics_conn:
            if code == stat[0]:
                mean = stat[1]
                std = stat[2]
                break
        # Probability of path type calculated through normal distr.
        probability_conn *= norm.pdf((path_num - mean) / std)
        # Check if the type is found in the disconnected statistics.
        # If not found, we assume mean = 0, and stdDev = 1.
        mean = 0
        std = 1
        for stat in statistics_disc:
            if code == stat[0]:
                mean = stat[1]
                std = stat[2]
                break
        # Probability of path type calculated through normal distr.
        probability_disc *= norm.pdf((path_num - stat[1]) / stat[2])
    return probability_conn, probability_disc
