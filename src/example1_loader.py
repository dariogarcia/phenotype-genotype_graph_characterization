"""
This file shows example on how to:
-Load the ontology data into data structures
-Load the data structures into a graph
-Find all paths between a pair of vertices
-Find all alternative paths between all pairs of phenotype-genotype connected through a gene
"""

from ontology_parser import load_data
from graph_methods import build_graph, find_all_paths
from extract_paths import get_connected_phenotype_genotype_alternative_paths
from extract_paths import get_disconnected_phenotype_genotype_paths

#OPTION 1: File paths obtained from arguments
#import sys
#load_data(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

#OPTION 2: File paths obtained from hardcoded constants
phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links = load_data('../data_files/hp.obo', '../data_files/go.obo', '../data_files/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt', '../data_files/goa_human.gaf')

#Check the sizes of the structures loaded
print 'Found',len(phenotypes),'phenotypes'
print 'Found',len(genotypes),'genotypes'
print 'Found',len(genes),'genes'
print 'Found',len(ph_ph_links,),'phenotype_phenotype links'
print 'Found',len(go_go_links,),'genotype_genotype links'
print 'Found',len(ph_gn_links,),'phenotype_gene links'
print 'Found',len(go_gn_links,),'genotype_gene links'

#Build a graph with the data loaded
#G = build_graph(phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links, undirected=True)

#Find all paths between a pair of vertices
#find_all_paths(G,G.vs.find(phenotypes[0]).index,G.vs.find(genotypes[0]).index,maxlen=4)

#Compute all alternative paths for connected phenotype-genotype pairs
#print 'Computing alternative paths among connected pairs'
#get_connected_phenotype_genotype_alternative_paths(phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links, continuing=True)

#Compute all paths for non-connected phenotype-genotype pairs
print 'Computing alternative paths among disconnected pairs'
get_disconnected_phenotype_genotype_paths(phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links, continuing=False)

