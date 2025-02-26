from src import dataset
from src import metrics
from src import plot

adata = dataset.mibitof()
G = dataset.spatial_graph(adata)

closeness = metrics.closeness(G)
pagerank = metrics.pagerank(G)

pos = adata.obsm['spatial']
edges = plot.edges(adata)
plot.Graph.plot(G, pos, edges, closeness, 'Mibitof with closeness centrality', './figures/mibitof/spatial_graph.svg')

hull = plot.convex_hull(pos)

plot.Quantification.plot(pos, closeness, 'Closeness centrality', './figures/mibitof/quantification/closeness.svg')
plot.Quantification.plot(pos, pagerank, 'Pagerank centrality', './figures/mibitof/quantification/pagerank.svg')
