# sq.pl.spatial_scatter(
#     adata,
#     shape=None,
#     # outline=True,
#     connectivity_key="spatial_connectivities",
#     library_id=["point8"],
#     seg_cell_id="cell_id",
#     color="Cluster",
#     library_key="library_id",
#     title=["point8"],
#     edges_color="white",
#     edges_width=0.1,
#     save='./mibitof/Spatial_scatter',
# )

# sq.pl.spatial_segment(
#     adata,
#     library_id=["point8"],
#     seg_cell_id="cell_id",
#     color="Cluster",
#     library_key="library_id",
#     title=["point8"],
#     save='./mibitof/Spacial_segment',
# )

import sys
import math
import numpy as np
import networkx as nx

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection

import squidpy as sq

# mibitof (non-grid dataset)
adata = sq.datasets.mibitof()

# planar graph of spatial neighbors
planar_spatial_connectivities, planar_spatial_distances = sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic", copy=True)

# connectivity matrix from spatial coordinates
sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
sq.gr.centrality_scores(adata, "Cluster", show_progress_bar=False)

# print(adata.uns['cluster_centrality_scores'])
# print(adata.obs['cluster']) # centrality scores are collected and merged for each "cluster" -> I want them for each node individually instead

adj = planar_spatial_connectivities.nonzero()
weight_adj = planar_spatial_distances.nonzero()

def cross_z(a, b):
    return a[0] * b[1] - a[1] * b[0]

def vec(start, end):
    return [end[0]-start[0],end[1]-start[1]]

def eql(point, other):
    if point is None or other is None:
        return False
    return point[0] == other[0] and point[1] == other[1]

def vec_len(start, end):
    vector = vec(start, end)
    return math.sqrt(vector[0]**2 + vector[1]**2)

spatial = adata.obsm['spatial']
spatial_connectivities = adata.obsp['spatial_connectivities']
x = spatial[:, 0]
y = spatial[:, 1]

min_x = min(x)
min_y = min(y)
max_x = max(x)
max_y = max(y)

edges = []
# -> map each index to a point in the plot
for start in range(0, len(x)):
    for end in range(0, len(y)):
        if spatial_connectivities[start, end] != 0:
            edges.append(LineCollection([np.column_stack([[x[start], x[end]],[y[start], y[end]]])], colors=['k'], linewidths=0.1))

# gift wrapping -> convex hull
convex_hull = []
hull_point = np.array([min_x,spatial[np.argmin(x)][1]])
convex_hull.append(hull_point)
endpoint = None
i = 0
while not eql(endpoint, convex_hull[0]):
    endpoint = spatial[0]
    for point in spatial:
        if eql(endpoint, hull_point) or cross_z(vec(convex_hull[i], endpoint), vec(convex_hull[i], point)) < 0:
            endpoint = point
    i += 1
    hull_point = endpoint
    convex_hull.append(hull_point)
# gift wrapping

nx_graph = nx.from_scipy_sparse_array(planar_spatial_distances)

# BEGIN TESTING
def draw(G, pos, measures, measure_name):
    
    nodes = nx.draw_networkx_nodes(G, pos, node_size=5, cmap=plt.cm.plasma, 
                                   node_color=list(measures.values()),
                                   nodelist=measures.keys())
    nodes.set_norm(mcolors.SymLogNorm(linthresh=0.01, linscale=1, base=10))
    # labels = nx.draw_networkx_labels(G, pos)
    edges = nx.draw_networkx_edges(G, pos)

    plt.title(measure_name)
    plt.colorbar(nodes)
    plt.axis('off')
    plt.show()

# draw(nx_graph, spatial, nx.closeness_centrality(nx_graph, distance='weight'), 'Closeness Centrality')
# END TESTING

# closeness
closeness_nodes = nx.closeness_centrality(nx_graph, distance='weight')
c = np.fromiter(closeness_nodes.values(), dtype=float)
closeness_max = max(closeness_nodes.values())
closeness_min = min(closeness_nodes.values())
# closeness

# pagerank
pagerank_nodes = nx.pagerank(nx_graph)
# pagerank

fig, ax = plt.subplots()
# ch = LineCollection([convex_hull], colors=['g'], linewidths=0.5)
# ax.add_collection(ch)
for edge in edges:
    ax.add_collection(edge)
sc = ax.scatter(x, y, s=1, c=c, cmap=plt.cm.plasma, vmin=closeness_min, vmax=closeness_max) # map closeness values as color mapping on the verticies
# ax.add_patch(patches.Rectangle((min_x, min_y), max_x - min_x, max_y - min_y, linewidth=0.2, edgecolor='r', facecolor="none"))
fig.colorbar(sc, ax=ax)
fig.savefig('./figures/mibitof/spatial_coordinates.svg', format='svg')

# min distance of each point in the spatial graph to the convex hull
distances = []
keys = iter(closeness_nodes.keys())
for point in spatial:
    min_distance = math.inf
    key = next(keys)
    for edge in convex_hull:
        distance = vec_len(point, edge)
        if distance < min_distance:
            min_distance = distance
    distances.append([key, min_distance, closeness_nodes[key], pagerank_nodes[key]])
# sort by distance
distances.sort(key=lambda entry: entry[1])
distances = np.array(distances)
# print(distances)

fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)
# Closeness-distance relationship
ax1.set_title('Closeness centrality')
ax1.set_xlabel('Distance to Bounding-Box')
ax1.set_ylabel('Centrality')
ax1.scatter(distances[:, 1], distances[:, 2], color='b', s=0.2)
# Pagerank-distance relationship
ax2.set_title('Pagerank centrality')
ax2.set_xlabel('Distance to Bounding-Box')
ax2.set_ylabel('Centrality')
ax2.scatter(distances[:, 1], distances[:, 3], color='c', s=0.2)
fig.savefig('./figures/mibitof/spatial_distances_to_hull.svg', format='svg')
