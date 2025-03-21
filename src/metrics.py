import networkx as nx

def closeness(G):
    """
    Calculate closeness metric using networkx.
    @param G [Graph] networkx Graph to calculate metric for.
    @return [Dict]
    """
    nodes = nx.closeness_centrality(G, distance='weight')
    return nodes

def pagerank(G, alpha=0.85):
    """
    Calculate pagerank metric using networkx.
    @param G [Graph] networkx Graph to calculate metric for.
    @param alpha [float] Value to use for pagerank algorithm. Value should be in [0,1]. *Default*: 0.85
    @return [Dict]
    """
    nodes = nx.pagerank(G, alpha=alpha)
    return nodes

# TODO: are other metrics also useful to have?
# - betweenness?
