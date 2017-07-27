"""
This file is similar to example 3, it shows an example on how to:
-Load the ontology data into a graph (as in example 1)
-Load statistics concerning different path types linking (dis)connected pairs
-Generate various disconnected pairs that were not part of the training set
 used to generate the statistics
-Use the estimator on each pair, and keep a record of the top 50 pairs with
 highest estimated probability of being connected
-Save the top 50 list (each pair and its estimated probability) using pickle
"""

from ontology_parser import *
from graph_methods import *
from extract_paths import *
from probability_estimator import *
import random
import signal
import sys

class GracefulKiller:
    """
    Safety measure to ensure that, even if the system kills the program,
    the generated data is not lost.
    Source: https://stackoverflow.com/questions/18499497/how-to-process-sigterm-signal-gracefully
    """
    kill_now = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self,signum, frame):
        self.kill_now = True


def find_top_50_pairs(phenotypes_ids, genotypes_ids, genes_ids,
                      phenotypes_links, genotypes_links,
                      phenotypes_genes_links, genotypes_genes_links,
                      total_results_path_conn='../results/total_results.pkl',
                      list_elems_path_conn='../results/list_elems.pkl',
                      type_index_path_conn='../results/type_index.pkl',
                      total_results_path_disc='../results/total_results_disc.pkl',
                      list_elems_path_disc='../results/list_elems_disc.pkl',
                      type_index_path_disc='../results/type_index_disc.pkl'):
    """
    Selects random phenotype-genotype pairs that are disconnected, until the
    program is killed or stopped with Ctrl+C, and stores in a list the top
    50 pairs with a highest estimated probability of being connected.
    This is useful to detect pairs which may be in fact connected, but are
    reported as disconnected because the gene that directly links them has
    not yet been found.

    Args:
        -phenotypes_ids: List of phenotypes
            +Type: list[str]
        -genotypes_ids: List of genotypes
            +Type: list[str]
        -genes_ids: List of genes
            +Type: list[str]
        -phenotypes_links: List of phenotype-phenotype links
            +Type: list[(str,str)]
        -genotypes_links: List of genotype-genotype links
            +Type: list[(str,str)]
        -phenotypes_genes_links: List of phenotype-genes links
            +Type: list[(str,str)]
        -genotypes_genes_links: List of genotype-genes links
            +Type: list[(str,str)]
        -total_results_path_conn: path to partial total_results for connected pairs
            +Type: str
        -list_elems_path_conn: path to partial list_elems for connected pairs
            +Type: str
        -type_index_path_conn: path to partial type_index for connected pairs
            +Type: str
        -total_results_path_disc: path to partial total_results for disconnected pairs
            +Type: str
        -list_elems_path_disc: path to partial list_elems for disconnected pairs
            +Type: str
        -type_index_path_disc: path to partial type_index for disconnected pairs
            +Type: str
    Returns:
        None. Persists data.
    """
    list_elems = pickle.load(open(list_elems_path_disc, 'rb'))
    graph = build_graph(phenotypes_ids, genotypes_ids, genes_ids, phenotypes_links, genotypes_links, phenotypes_genes_links, genotypes_genes_links)
    statistics_conn = get_stats_from_alternative_paths(total_results_path_conn, type_index_path_conn, list_elems_path_conn)
    statistics_disc = get_stats_from_alternative_paths(total_results_path_disc, type_index_path_disc, list_elems_path_disc)
    top_50_list = []
    prob_conn = prob_disc = 0
    loop_killer = GracefulKiller()

    try:
        # Loop until sufficient number of pairs have been tested
        while True:
            # Pick random phenotype
            p_id = random.choice(phenotypes_ids)
            # Get the list of linked genes
            p_genes = list(set([i[1] for i in phenotypes_genes_links if i[0] == p_id]))
            # And the list of linked genotypes
            p_linked_genotypes = list(set([i[0] for i in genotypes_genes_links if i[1] in p_genes]))
            # From all genotypes, take away the linked genotypes to get the disconnected genotypes
            p_genotypes = list(set(genotypes_ids).difference(p_linked_genotypes))
            # If the list of genotypes is empty, we start this itteration over
            if len(p_genotypes) == 0: continue

            # Pick a random genotype that is disconnected to the phenotype
            g_id = random.choice(p_genotypes)
            # If the selected pair was already in list_elems, we start this itteration over
            if (p_id, g_id) in list_elems: continue

            # Estimate probabilities the phenotype-genotype pair
            prob_conn, prob_disc = estimate_probabilities(graph, [], (p_id, g_id), statistics_conn, statistics_disc)

            # Append the results to the probabilities_table
            probabilities_table.append(((p_id, g_id), prob_conn, prob_disc, prob_diff))

            # Exit the loop if the program is killed
            if loop_killer.kill_now:
                break
        raise KeyboardInterrupt

    # Code executed at the end of the function, or if the program is interrupted/killed
    except KeyboardInterrupt:
        # Sort probabilities_table by probability difference
        probabilities_table = sorted(probabilities_table, key=lambda x: x[3])
        pickle.dump(probabilities_table, open("../results/estimated_probs_disc.pkl", "wb"))

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------


# Load the ontology data into data structures
print "Loading data into structures..."
phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links = load_data('../data_files/hp.obo', '../data_files/go.obo', '../data_files/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt', '../data_files/goa_human.gaf')
print "DONE!\n"

# Find the top 50 most likely disconnected pairs to be in fact connected
print "Finding the top 50 disconnected pairs..."
find_top_50_pairs(phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links)
print "DONE!\n"
