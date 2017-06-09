"""
This file contains methods related with the graph representation
and processing of the data loaded from the phenotype and 
genotype ontologies.
"""


def get_phenotype_genotype_path_code(G,path,p_id,common_genes,genotype_id):
    """
    Given a path in the graph, compute its code defined by a series of characters.
    P: Source phenotype
    p: other phenotype
    N: Common gene (link P-N has been removed from the graph)
    n: other gene
    G: Target genotype
    g: other genotype

    Args:
        -G: Graph where the path was computed.
            +Type: igraph graph
        -path: Path to be codified
            +Type: list[int]
        -p_id: Phenotype source of the path
            +Type: str
        -common_genes: List of genes linked with source phenotype removed from G
            +Type: list[str]
        -genotype_id: Genotype target of the path
            +Type: str
    Returns:
        -code: Code corresponding to the input path
            +Type: str
    """
    #Fill the code character by character
    code = ''
    for p in path:
        #If its a phenotype
        if G.vs[p]['name'][:3]=='HP:':
            #If its the source phenotype
            if G.vs[p]['name']==p_id:
                code = code + 'P' 
            else:
                code = code + 'p' 
        #If its a genotype
        elif G.vs[p]['name'][:3]=='GO:':
            #If its the target genotype
            if G.vs[p]['name']==genotype_id:
                code = code + 'G' 
            else:
                code = code + 'g' 
        #If its a gene
        else: 
            #If its a common gene
            if G.vs[p]['name'] in common_genes:
                code = code + 'N' 
            else:
                code = code + 'n' 
    return code


def build_graph(phenotypes_ids, genotypes_ids, genes_ids, phenotypes_links, genotypes_links, phenotypes_genes_links, genotypes_genes_links, undirected = True):
    """
    Generate a graph from a set of phenotypes, genotypes and genes, and links between them

    Args:
        -phenotypes_ids: List of phenotypes to be used as vertices
            +Type: list[str]
        -genotypes_ids: List of genotypes to be used as vertices
            +Type: list[str]
        -genes_ids: List of genes to be used as vertices
            +Type: list[str]
        -phenotypes_links: List of phenotype-phenotype links to be used as edges
            +Type: list[(str,str)]
        -genotypes_links: List of genotype-genotype links to be used as edges
            +Type: list[(str,str)]
        -phenotypes_genes_links: List of phenotype-genes links to be used as edges
            +Type: list[(str,str)]
        -genotypes_genes_links: List of genotype-genes links to be used as edges
            +Type: list[(str,str)]
    Returns:
        -G: Built graph
            +Type: igraph graph
    """
    import igraph as ig
    #If the graph is undirected
    if undirected:
        G = ig.Graph(directed=False)
        #Load vertices
        G.add_vertices(phenotypes_ids)
        G.add_vertices(genes_ids)
        G.add_vertices(genotypes_ids)
        #Load edges
        G.add_edges(phenotypes_links)
        G.add_edges(genotypes_links)
        G.add_edges(phenotypes_genes_links)
        G.add_edges(genotypes_genes_links)
    #If the graph is directed
    else: 
        G = ig.Graph(directed=True)
        #Load vertices
        G.add_vertices(phenotypes_ids)
        G.add_vertices(genes_ids)
        G.add_vertices(genotypes_ids)
        #Load edges
        G.add_edges(phenotypes_links)
        G.add_edges(genotypes_links)
        G.add_edges(phenotypes_genes_links)
        G.add_edges(genotypes_genes_links)
    #Return the graph
    return G


def find_all_paths(graph, start, end, mode = 'ALL', maxlen = None):
    """
    Find all paths of a maximum length between two 
    vertices of a graph

    Args:
        -graph: Graph where paths are to be found
            +Type: igraph graph
        -start: Source vertex index of the path
            +Type: int
        -end: Target vertex index of the path
            +Type: int
        -mode: Explore OUT/IN/ALL edges for building the paths 
            +Type: str
        -maxlen: Maximum path length to be xplored
            +Type: int
    Returns:
        -all_paths: list of paths from start to end
            +Type: list[list[int]]
    Source:
        https://stackoverflow.com/questions/29314795/python-igraph-get-all-possible-paths-in-a-directed-graph
        https://stackoverflow.com/questions/2606018/path-between-two-nodes
    """
    #Recursive method
    def find_all_paths_aux(adjlist, start, end, path, maxlen = None):
        #Add current vertex to the path
        path = path + [start]
        #Check if path is ended
        if start == end:
            return [path]
        paths = []
        #If length is under the maximum
        if maxlen is None or len(path) <= maxlen:
            #Make recursive call for every neighbour not yet in the path
            for node in adjlist[start] - set(path):
                paths.extend(find_all_paths_aux(adjlist, node, end, path, maxlen))
        return paths
    #Build the adjacency list of the graph
    adjlist = [set(graph.neighbors(node, mode = mode)) \
        for node in xrange(graph.vcount())]
    #Initiate recursive call
    all_paths = find_all_paths_aux(adjlist, start, end, [], maxlen)
    #Returns list of all paths
    return all_paths
