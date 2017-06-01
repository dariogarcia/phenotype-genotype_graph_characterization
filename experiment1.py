from __future__ import division
from collections import Counter
from multiprocessing import Pool, cpu_count
import cPickle as pickle
#Load data
from graph_generator_nx import *

def generate_path_code(G,path,p_id,common_genes,genotype_id):
    code = ''
    for p in path:
        if G.vs[p]['name'][:3]=='HP:':
            if G.vs[p]['name']==p_id:
                code = code + 'P'
            else:
                code = code + 'p'
        elif G.vs[p]['name'][:3]=='GO:':
            if G.vs[p]['name']==genotype_id:
                code = code + 'G'
            else:
                code = code + 'g'
        else: 
            if G.vs[p]['name'] in common_genes:
                code = code + 'N'
            else:
                code = code + 'n'
    return code
   
def compute_phenotype_paths(p_id):
    #build graph
    G_base = create_igraph_graph(phenotypes_ids, genotypes_ids, all_genes_ids, phenotypes_links, genotypes_links, genes_phenotypes_links, genes_genotypes_links, undirected=True)
    #prepare structures storing the results
    connected_genotypes_data = Counter()
    counter_connected = 0
    unconnected_genotypes_data = Counter()
    counter_unconnected = 0
    #precompute the set of genes connected with the current phenotype
    phenotype_genes = set([i[1] for i in genes_phenotypes_links if i[0]==p_id])
    #For each genotype
    for g_id in genotypes_ids:
        #find the shared genes
        genotype_genes = [i[0] for i in genes_genotypes_links if i[1]==g_id]
        common_gens = phenotype_genes.intersection(set(genotype_genes))
        #compute all the links to be removed and remove them
        to_be_removed_links = []
        for gen_id in common_gens:
            to_be_removed_links.append((p_id,gen_id))
        genes_phenotypes_links_pruned = [i for i in genes_phenotypes_links if i not in to_be_removed_links]
        #build graph if necessary
        if len(common_gens)>0:
            G = create_igraph_graph(phenotypes_ids, genotypes_ids, all_genes_ids, phenotypes_links, genotypes_links, genes_phenotypes_links_pruned, genes_genotypes_links, undirected=True)
        else:
            G = G_base
        paths_codes = Counter()
        paths = find_all_paths(G,G.vs.find(p_id).index,G.vs.find(g_id).index,maxlen=4)
        #Explore the statistics of each path
        for current_path in paths:
            current_code  = generate_path_code(G,current_path,p_id,common_gens,g_id)
            if current_code in paths_codes:
                paths_codes[current_code]+=1
            else:
                paths_codes[current_code]=1
        #Update the dictionaries with the paths info
        if len(common_gens)>0:
            counter_connected+=1
            connected_genotypes_data+=paths_codes
        else:
            counter_unconnected+=1
            unconnected_genotypes_data+=paths_codes
    #Normalize values
    for key, value in connected_genotypes_data.items():
            connected_genotypes_data[key] = value / counter_connected
    for key, value in unconnected_genotypes_data.items():
            unconnected_genotypes_data[key] = value / counter_unconnected
    pickle.dump(connected_genotypes_data, open( "connected_genotypes_data_"+str(p_id)+".pkl", "wb" ) )
    pickle.dump(unconnected_genotypes_data, open( "unconnected_genotypes_data_"+str(p_id)+".pkl", "wb" ) )
    return connected_genotypes_data,unconnected_genotypes_data



phenotypes_ids, genotypes_ids, all_genes_ids, phenotypes_links, genotypes_links, genes_phenotypes_links, genes_genotypes_links = load_consistent_structures()

total_connected_genotypes_data = Counter()
total_unconnected_genotypes_data = Counter()
pool = Pool(cpu_count())
list_genotypes = pool.map(compute_phenotype_paths, phenotypes_ids)
#Aggregate all values
for x in list_genotypes:
    total_connected_genotypes_data+=x[0]
    total_unconnected_genotypes_data+=x[1]
#Normalize
for key, value in total_connected_genotypes_data.items():
    total_connected_genotypes_data[key] = value / len(phenotypes_ids)
for key, value in total_unconnected_genotypes_data.items():
    total_unconnected_genotypes_data[key] = value / len(phenotypes_ids)
#Store
pickle.dump( total_connected_genotypes_data, open( "total_connected_genotypes_data.pkl", "wb" ) )
pickle.dump( total_unconnected_genotypes_data, open( "total_unconnected_genotypes_data.pkl", "wb" ) )
print 'FINAL TOTAL C',total_connected_genotypes_data
print 'FINAL TOTAL UNC',total_unconnected_genotypes_data


##For each gene connected
#for gene_id in linked_genes:
#    #Get list of genotypes connected to gene
#    linked_genotypes = [i[1] for i in genes_genotypes_links if i[0]==gene_id]
#    #Remove link
#    genes_phenotypes_links_pruned = [i for i in genes_phenotypes_links if i != (p_id,gene_id)]
#    #Build the graph
#    G = create_igraph_graph(phenotypes_ids, genotypes_ids, all_genes_ids, phenotypes_links, genotypes_links, genes_phenotypes_links_pruned, genes_genotypes_links, undirected=True)
#    #For each genotype previously linked through the removed gene, explore the paths
#    for genotype_id in linked_genotypes:
#        paths_codes = {}
#        paths = find_all_paths(G,G.vs.find(p_id).index,G.vs.find(genotype_id).index,maxlen=5)
#        #Explore the statistics of each path
#        for current_path in paths:
#            current_code  = generate_path_code(G,current_path,p_id,gene_id,genotype_id)
#            if current_code in paths_codes:
#                paths_codes[current_code]+=1
#            else:
#                paths_codes[current_code]=1
#            if 'PG' in current_code:
#                print 'Impossible path',current_path
#        print p_id, genotype_id, paths_codes





