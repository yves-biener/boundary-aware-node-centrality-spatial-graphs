from src import dataset
from src import metrics
from src import plot

adata = dataset.merfish()
G = dataset.spatial_graph(adata)

closeness = metrics.closeness(G)
pagerank = metrics.pagerank(G)

pos = adata.obsm['spatial']
edges = plot.edges(adata)
plot.Graph.plot(G, pos, edges, closeness, 'Merfish with closeness centrality', './figures/merfish/spatial_graph.svg')

hull = plot.convex_hull(pos)

plot.Quantification.plot(pos, closeness, 'Closeness centrality', './figures/merfish/quantification/closeness.svg')
plot.Quantification.plot(pos, pagerank, 'Pagerank centrality', './figures/merfish/quantification/pagerank.svg')
