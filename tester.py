from ontology_parser import *

#import sys
#load_data(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
phenotypes, genotypes, genes, ph_ph_links, go_go_links, ph_gn_links, go_gn_links = load_data('../data_files/hp.obo', '../data_files/go.obo', '../data_files/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt', '../data_files/goa_human.gaf')

print 'Found',len(phenotypes),'phenotypes'
print 'Found',len(genotypes),'genotypes'
print 'Found',len(genes),'genes'
print 'Found',len(ph_ph_links,),'phenotype_phenotype links'
print 'Found',len(go_go_links,),'genotype_genotype links'
print 'Found',len(ph_gn_links,),'phenotype_gene links'
print 'Found',len(go_gn_links,),'genotype_gene links'
