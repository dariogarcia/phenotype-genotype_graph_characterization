

#Read phenotype Ids from Phenotype Ontology
#Return list of unique phenotypes Ids
def load_phenotypes():
    import os
    path = '../data_files/'
    phenotypes_ids = []
    with open(os.path.join(path,'hp.obo')) as f:
        for line in f:
            #Found a new term
            if line=='[Term]\n':
                term_id = next(f).split()[1]
                if 'HP' not in term_id:
                    raise Exception('Wrong phenotype id')
                phenotypes_ids.append(term_id.splitlines()[0])
    return list(set(phenotypes_ids))


#Read phenotype-phenotype is_a links from Phenotype Ontology
#Returns list of tuples (source_phenotypeId,target_phenotypeId)
def load_phenotype_links():
    import os
    path = '../data_files/'
    phenotypes_links = []
    with open(os.path.join(path,'hp.obo')) as f:
        source_phenotype = ''
        target_phenotype = ''
        for line in f:
            #Found a new term
            if line.startswith('id: '):
                source_phenotype = line.split()[1].splitlines()[0]
                if 'HP' not in source_phenotype:
                    raise Exception('Wrong source phenotype id')
                continue
            #Make sure we have a source at this point
            if source_phenotype == ' ':
                raise Exception('Undefined source phenotype id')
            #Find the target term
            if line.startswith('is_a: '):
                target_phenotype = line.split()[1].splitlines()[0]
                if 'HP' not in target_phenotype:
                    raise Exception('Wrong target phenotype id')
                phenotypes_links.append((source_phenotype,target_phenotype))
                continue
    return list(set(phenotypes_links))

    
#Read genes Ids from phenotypes-gene and genotype-gene links
#Returns list of unique genes Ids
def load_all_genes():
    import os
    path = '../data_files/'
    #Read genes from phenotype-gene links
    genes_ids = set()
    with open(os.path.join(path,'ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt')) as f:
        #Skip the first line of syntax
        next(f)
        for line in f:
            gene_id = line.split()[1]
            if gene_id not in genes_ids:
                genes_ids.add(gene_id)
    #Read genes from genotype-gene links
    genes_ids2 = set()
    with open('../data_files/goa_human.gaf') as f:
        for line in f:
            #Skip comment lines
            if line[0] == '!':
                continue
            gene_id = line.split()[2]
            genes_ids2.add(gene_id)
    return list(genes_ids.union(genes_ids2))

    
#Get the list of genes Ids which are shared by GO and PHO links
#Returns list of gene Ids
def load_genes_from_links():
    genes_phenotypes_links = load_genes_phenotypes_links()
    genes_phenotypes = set([i[1] for i in genes_phenotypes_links])
    genes_genotypes_links = load_genes_genotypes_links()
    genes_genotypes = set([i[0] for i in genes_genotypes_links])
    return list(genes_genotypes.intersection(genes_phenotypes))


#Read links between genes and phenotypes
#Returns tuples of (geneId,phenotypeId)
def load_genes_phenotypes_links():
    import os
    path = '../data_files/'
    genes_phenotypes_links = []
    with open(os.path.join(path,'ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt')) as f:
        #Skip the first line of syntax
        next(f)
        for line in f:
            source_gene = line.split()[1]
            target_phenotype = line.split()[-1:][0]
            if 'HP' not in target_phenotype:
                raise Exception('Wrong target phenotype id for link')
            #genes_phenotypes_links.append((source_gene,target_phenotype))
            genes_phenotypes_links.append((target_phenotype,source_gene))
    return list(set(genes_phenotypes_links))

#Read genotypes Ids from GO
#Returns list of unique genotypes Ids
def load_genotypes():
    import os
    path = '../data_files/'
    genotypes_ids = []
    with open(os.path.join(path,'go.obo')) as f:
        for line in f:
            #Found a new term
            if line=='[Term]\n':
                genotypes_ids.append(next(f).split()[1].splitlines()[0])
    return list(set(genotypes_ids))
   

#Read geontype-genotype is_a links from the GO
#Returns list of tuples (source_genotypeId,target_genotypeId)
def load_genotype_links():
    import os
    path = '../data_files/'
    genotypes_links = []
    with open(os.path.join(path,'go.obo')) as f:
        source_genotype = ''
        target_genotype = ''
        for line in f:
            #Stop reading when all terms have been read
            if line.startswith('[Typedef]'):
                break
            #Found a new term
            if line.startswith('id: '):
                source_genotype = line.split()[1].splitlines()[0]
                if 'GO' not in source_genotype:
                    raise Exception('Wrong source genotype id')
                continue
            #Make sure we have a source at this point
            if source_genotype == ' ':
                raise Exception('undefined source genotype id')
            #Find the target term
            if line.startswith('is_a: '):
                target_genotype = line.split()[1].splitlines()[0]
                if 'GO' not in target_genotype:
                    raise Exception('Wrong target genotype id')
                genotypes_links.append((source_genotype,target_genotype))
                continue
    return list(set(genotypes_links))
    
#Read links between genes and genotypes
#Return list of tuples (geneId,genotypeId)
def load_genes_genotypes_links():
    import os
    path = '../data_files/'
    genes_genotypes_links = []
    with open(os.path.join(path,'goa_human.gaf')) as f:
        for line in f:
            #Skip comment lines
            if line[0] == '!':
                continue
            #In some cases the split needs an offset
            offset = 0
            if line.split()[6] in ['EXP','IDA','IPI','IMP','IGI','IEP','TAS','ND','NAS','ISS','IKR','IEA','IC','IBA']:
                offset = 1
            #Avoid those links which are not experimental results
            if line.split()[5+offset] not in ['EXP','IDA','IPI','IMP','IGI','IEP']:
                continue
            source_gene = line.split()[2+offset]
            target_genotype = line.split()[3+offset]
            if 'GO:' not in target_genotype:
                raise Exception('Wrong target genotype id in link')
            genes_genotypes_links.append((source_gene,target_genotype))
    return list(set(genes_genotypes_links))
   
#Loads and returns phenotypes, genotypes, genes, and their links
#Only genes shared by both ontologies, and links with those genes, are returned
def load_consistent_structures():
    g2 = load_genes_from_links()
    #Avoid links which relate genes not linked with both ontologies
    gp = load_genes_phenotypes_links()
    clean_gp_links = [l for l in gp if l[1] in g2]
    ggt = load_genes_genotypes_links()
    clean_ggt_links = [l for l in ggt if l[0] in g2]
    return load_phenotypes(), load_genotypes(), g2, load_phenotype_links(), load_genotype_links(), clean_gp_links, clean_ggt_links 

#phenotypes_ids, genotypes_ids, all_genes_ids, phenotypes_links, genotypes_links, genes_phenotypes_links, genes_genotypes_links = load_consistent_structures()

def create_igraph_graph(phenotypes_ids, genotypes_ids, all_genes_ids, phenotypes_links, genotypes_links, genes_phenotypes_links, genes_genotypes_links,Pgd = 0, PPd = 0, Ggd = 0, GGd = 0, undirected = False):
    #Parse ontologies
    #Initialize graph
    import igraph as ig
    if undirected:
        G = ig.Graph(directed=False)
        #Load vertices
        G.add_vertices(phenotypes_ids)
        G.add_vertices(all_genes_ids)
        G.add_vertices(genotypes_ids)
        #Load edges
        G.add_edges(phenotypes_links)
        G.add_edges(genotypes_links)
        G.add_edges(genes_phenotypes_links)
        G.add_edges(genes_genotypes_links)
    else: 
        G = ig.Graph(directed=True)
        #Load vertices
        G.add_vertices(phenotypes_ids)
        G.add_vertices(all_genes_ids)
        G.add_vertices(genotypes_ids)
        #Load edges
        G.add_edges(phenotypes_links)
        G.add_edges(genotypes_links)
        G.add_edges(genes_phenotypes_links)
        G.add_edges(genes_genotypes_links)
    return G

def find_all_paths(graph, start, end, mode = 'OUT', maxlen = None):
    def find_all_paths_aux(adjlist, start, end, path, maxlen = None):
        path = path + [start]
        if start == end:
            return [path]
        paths = []
        if maxlen is None or len(path) <= maxlen:
            for node in adjlist[start] - set(path):
                paths.extend(find_all_paths_aux(adjlist, node, end, path, maxlen))
        return paths
    adjlist = [set(graph.neighbors(node, mode = mode)) \
        for node in xrange(graph.vcount())]
    all_paths = []
    start = start if type(start) is list else [start]
    end = end if type(end) is list else [end]
    for s in start:
        for e in end:
            all_paths.extend(find_all_paths_aux(adjlist, s, e, [], maxlen))
    return all_paths


#Creates a networkx graph
def create_networkx_graph(phenotypes_ids, genotypes_ids, all_genes_ids, phenotypes_links, genotypes_links, genes_phenotypes_links, genes_genotypes_links):
    #Parse ontologies
    #Initialize graph
    import networkx as nx
    G=nx.Graph()
    #Load vertices
    G.add_nodes_from(phenotypes_ids, type='phenotype')
    G.add_nodes_from(all_genes_ids, type='gene')
    G.add_nodes_from(genotypes_ids, type='genotype')
    #Load edges
    G.add_edges_from(phenotypes_links)
    G.add_edges_from(genotypes_links)
    G.add_edges_from(genes_phenotypes_links)
    G.add_edges_from(genes_genotypes_links)
    return G
