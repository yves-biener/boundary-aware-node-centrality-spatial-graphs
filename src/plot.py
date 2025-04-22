import math

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from matplotlib.collections import LineCollection

class Vector:
    """
    Small Vector class used for calculations of distances
    """
    def cross_z(a, b):
        return a[0] * b[1] - a[1] * b[0]

    def vec(start, end):
        """
        Create a new vector from `start` to `end`
        @param start [Array-2d]
        @param end [Array-2d]
        @return [Array-2d]
        """
        return [end[0]-start[0],end[1]-start[1]]

    def eql(point, other):
        """
        Check if `point` is equivalent to `other`.
        """
        if point is None or other is None:
            return False
        return point[0] == other[0] and point[1] == other[1]

    def vec_len(vector):
        """
        Calculate the length of a `vector`
        """
        return math.sqrt(vector[0]**2 + vector[1]**2)


def edges(adata):
    """
    Create an array of `LineCollection`'s which describe the edges of the spatial graph inside of `adata`.
    """
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
    """
    Calculate the convex hull using gift wrapping.
    @param pos [Array-2d]
    @return [Array-2d] ordered point array of hull points
    """
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
    """
    Graph class containing functions for plotting
    """
    def plot(G, pos, edges, measures, name, path):
        """
        Plot a given graph.
        @param G [Graph] networkx Graph
        @param pos [Array-2d] position values of points
        @param edges [Array] LineCollection's that should be used to render the edges. @see edges(pos)
        @param measures [Dict] Result of a centrality calculation of networkx on graph G that should be visualized
        @param name [String] Name of the measurement used, which will be integrated into the generated plot
        @param path [String] Path to store the generated plot as svg file
        """
        fig, ax = plt.subplots()
        x = pos[:, 0]
        y = pos[:, 1]
        c = np.fromiter(measures.values(), dtype=float)
        # convex hull -> Bounding-Box
        ch = LineCollection([convex_hull(pos)], colors=['r'], linewidths=1)
        ax.add_collection(ch)
        for edge in edges:
            ax.add_collection(edge)
        sc = ax.scatter(x, y, s=1, cmap=plt.cm.plasma, c=c) # map closeness values as color mapping on the verticies
        ax.set_title(name)
        fig.colorbar(sc, ax=ax)
        fig.savefig(path, format='svg')

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


def normalize_dict(d):
    max = np.max(list(d.values()))
    return {k: (v / max) for k, v in d.items()}


class Quantification:
    def data(pos, metric):
        """
        Create distance to metric relationship.
        @param pos [Array-2d] Positions of points
        @param metric [Dict] Result of a centrality calculation of networkx on graph G
        @return [Array-2d] relationship ordered (acending) by distance
        """
        quantification = []
        keys = iter(metric.keys())
        hull = convex_hull(pos)

        for point in pos:
            min_distance = math.inf
            key = next(keys)
            for edge in hull:
                vector = Vector.vec(point, edge)
                distance = Vector.vec_len(vector)
                if distance < min_distance:
                    min_distance = distance
            quantification.append([min_distance, metric[key]])

        # sort by distance
        quantification.sort(key=lambda entry: entry[0])
        return np.array(quantification)

    def opt_data(pos, metric, b_opt):
        """
        Optimal data relationship. Describe Points which are effected by the boundary.
        @param pos [Array-2d] Positions of points
        @param metric [Dict] Result of a centrality calculation of networkx on graph G
        @param b_opt [int] Cross point distance.
        @return [Dict] `metric` like dict (same keys) where True determines a node that is uneffected by the boundary, otherwise False.
        """
        keys = iter(metric.keys())
        boundary_effected = {}
        hull = convex_hull(pos)

        for point in pos:
            min_distance = math.inf
            key = next(keys)
            for edge in hull:
                vector = Vector.vec(point, edge)
                distance = Vector.vec_len(vector)
                if distance < min_distance:
                    min_distance = distance
            if min_distance <= b_opt:
                boundary_effected[key] = True
            else:
                boundary_effected[key] = False

        # sort by distance
        return boundary_effected

    def plot(data, d_curve, C_curve, metric_name, path):
        """
        Plot relationship data.
        @param data [Array-2d] see `data(pos, metric)`
        @param d_curve linear function of the left side of the intersection point
        @param C_curve constant function of the right side of the intersection point
        @param metric_name [String] Name of the metric to be used as a title for the plot
        @param path [String] Path to store the generated plot as svg file
        """
        fig, ax = plt.subplots()
        ax.set_title(metric_name)
        ax.set_xlabel('Distance to Bounding-Box')
        ax.set_ylabel('Centrality')
        ax.scatter(data[:, 0], data[:, 1], color='b', s=0.2)
        if d_curve is not None and C_curve is not None:
            ax.plot(d_curve, C_curve, color='r', linewidth=1)
        fig.savefig(path, format='svg')
