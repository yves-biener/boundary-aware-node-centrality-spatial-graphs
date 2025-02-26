import math

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from matplotlib.collections import LineCollection

class Vector:
    def cross_z(a, b):
        return a[0] * b[1] - a[1] * b[0]

    def vec(start, end):
        return [end[0]-start[0],end[1]-start[1]]

    def eql(point, other):
        if point is None or other is None:
            return False
        return point[0] == other[0] and point[1] == other[1]

    def vec_len(start, end):
        vector = Vector.vec(start, end)
        return math.sqrt(vector[0]**2 + vector[1]**2)


def edges(adata):
    spatial = adata.obsm['spatial']
    spatial_connectivities = adata.obsp['spatial_connectivities']
    x = spatial[:, 0]
    y = spatial[:, 1]

    e = []
    for start in range(0, len(x)):
        for end in range(0, len(y)):
            if spatial_connectivities[start, end] != 0:
                e.append(LineCollection([np.column_stack([[x[start], x[end]],[y[start], y[end]]])], colors=['k'], linewidths=0.1))
    return e


def convex_hull(pos):
    # gift wrapping algorithm
    x = pos[:, 0]
    y = pos[:, 1]

    min_x = min(x)
    min_y = min(y)
    max_x = max(x)
    max_y = max(y)

    convex_hull = []
    hull_point = np.array([min_x,pos[np.argmin(x)][1]])
    convex_hull.append(hull_point)
    endpoint = None
    i = 0
    while not Vector.eql(endpoint, convex_hull[0]):
        endpoint = pos[0]
        for point in pos:
            if Vector.eql(endpoint, hull_point) or Vector.cross_z(Vector.vec(convex_hull[i], endpoint), Vector.vec(convex_hull[i], point)) < 0:
                endpoint = point
        i += 1
        hull_point = endpoint
        convex_hull.append(hull_point)
    return convex_hull


class Graph:
    def plot(G, pos, edges, measures, name, path, **kwargs):
        fig, ax = plt.subplots()
        x = pos[:, 0]
        y = pos[:, 1]
        c = np.fromiter(measures.values(), dtype=float)
        # TODO: get optional hull to draw from optional **kwargs argument
        # ch = LineCollection([convex_hull], colors=['g'], linewidths=0.5)
        # ax.add_collection(ch)
        for edge in edges:
            ax.add_collection(edge)
        sc = ax.scatter(x, y, s=1, cmap=plt.cm.plasma, c=c) # map closeness values as color mapping on the verticies
        # TODO: draw bounding box?
        # ax.add_patch(patches.Rectangle((min_x, min_y), max_x - min_x, max_y - min_y, linewidth=0.2, edgecolor='r', facecolor="none"))
        fig.colorbar(sc, ax=ax)
        fig.savefig('name', format='svg')

    def show(G, pos, measures, name):
        nodes = nx.draw_networkx_nodes(G, pos, node_size=2, cmap=plt.cm.plasma, 
                                       node_color=list(measures.values()),
                                       nodelist=measures.keys())
        nodes.set_norm(mcolors.SymLogNorm(linthresh=0.01, linscale=1, base=10))
        edges = nx.draw_networkx_edges(G, pos)

        plt.title(name)
        plt.colorbar(nodes)
        plt.axis('off')
        plt.show()


class Quantification:
    def plot(pos, metric, metric_name, path):
        quantification = []

        hull = convex_hull(pos)
        keys = iter(metric.keys())

        for point in pos:
            min_distance = math.inf
            key = next(keys)
            for edge in hull:
                distance = Vector.vec_len(point, edge)
                if distance < min_distance:
                    min_distance = distance
            quantification.append([min_distance, metric[key]])

        # sort by distance
        quantification.sort(key=lambda entry: entry[0])
        quantification = np.array(quantification)

        fig, ax = plt.subplots()
        ax.set_title(metric_name)
        ax.set_xlabel('Distance to Bounding-Box')
        ax.set_ylabel('Centrality')
        ax.scatter(quantification[:, 0], quantification[:, 1], color='b', s=0.2)
        fig.savefig(path, format='svg')

