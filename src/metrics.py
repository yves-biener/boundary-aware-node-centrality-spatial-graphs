import networkx as nx

def closeness(G):
    nodes = nx.closeness_centrality(G, distance='weight')
    return nodes

def pagerank(G):
    nodes = nx.pagerank(G)
    return nodes

# TODO: are other metrics also useful to have?
