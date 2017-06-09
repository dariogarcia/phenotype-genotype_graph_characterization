from ontology_parser import *
from graph_methods import *

#Load the files from input arguments
#import sys
#load_data(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

#Or from hardcoded paths
phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links = load_data('../data_files/hp.obo', '../data_files/go.obo', '../data_files/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt', '../data_files/goa_human.gaf')

#Check sizes
#print 'Found',len(phenotypes),'phenotypes'
#print 'Found',len(genotypes),'genotypes'
#print 'Found',len(genes),'genes'
#print 'Found',len(ph_ph_links,),'phenotype_phenotype links'
#print 'Found',len(go_go_links,),'genotype_genotype links'
#print 'Found',len(ph_gn_links,),'phenotype_gene links'
#print 'Found',len(go_gn_links,),'genotype_gene links'

#Build a graph with it
G = build_graph(phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links, undirected=True)

#Find all paths between a pair of vertices
find_all_paths(G,G.vs.find(phenotypes[0]).index,G.vs.find(genotypes[0]).index,maxlen=4)

