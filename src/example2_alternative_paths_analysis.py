"""
This file shows an example on how to process the alternative paths
previously computed using "get_connected_phenotype_genotype_alternative_paths"
"""

# from graph_methods import * # This isn't used anymore
from extract_paths import read_and_analyze_alternative_paths
# from extract_paths import get_connected_phenotype_genotype_alternative_paths
# from extract_paths import get_disconnected_phenotype_genotype_paths

#This code assumes the alternative paths have already been computed (at least partially)
#using the following method (see example1_loader.py for details)
#get_connected_phenotype_genotype_alternative_paths(phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links)


#OPTION 1: File paths obtained from arguments
#import sys
#read_and_analyze_alternative_paths(sys.argv[1],sys.argv[2],sys.argv[3])

#OPTION 2: File paths hardcoded into constants
#print 'Connected pairs'
#read_and_analyze_alternative_paths('../results/total_results.pkl', '../results/type_index.pkl', '../results/list_elems.pkl')
#print 'Disconnected pairs'
#read_and_analyze_alternative_paths('../results/total_results_disc.pkl', '../results/type_index_disc.pkl', '../results/list_elems_disc.pkl')
print 'Connected pairs'
read_and_analyze_alternative_paths('../results/old_versions/total_results.pkl', '../results/old_versions/type_index.pkl', '../results/old_versions/list_elems.pkl')
print 'Disconnected pairs'
read_and_analyze_alternative_paths('../results/old_versions/total_results_disc.pkl', '../results/old_versions/type_index_disc.pkl', '../results/old_versions/list_elems_disc.pkl')
