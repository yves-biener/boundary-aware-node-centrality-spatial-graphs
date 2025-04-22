import math
import networkx as nx

from src import plot

def closeness(G):
    """
    Calculate closeness metric using networkx.
    @param G [Graph] networkx Graph to calculate metric for.
    @return [Dict] with normalized values
    """
    nodes = nx.closeness_centrality(G, distance='weight')
    return plot.normalize_dict(nodes)

def pagerank(G, alpha=0.85):
    """
    Calculate pagerank metric using networkx.
    @param G [Graph] networkx Graph to calculate metric for.
    @param alpha [float] Value to use for pagerank algorithm. Value should be in [0,1]. *Default*: 0.85
    @return [Dict] with normalized values
    """
    nodes = nx.pagerank(G, alpha=alpha)
    return plot.normalize_dict(nodes)

# TODO: are other metrics also useful to have?
# - betweenness?


class Correction:
    def correct(pos, metric, m_opt, c0_opt, b_opt):
        """
        Correct a given metric by the calculated model.
        @param pos [Array-2d] Positions of points
        @param metric [Dict] Result of a centrality calculation of networkx on graph G
        @param m_opt [Float] Model m value (slope of linear function)
        @param c0_opt [Float] Model c0 value (offset of linear function)
        @param b_opt [Float] Model b value (intersection point)
        @return [Dict] Corrected centrality values based on @param metric
        """
        corrected_metric = {}
        keys = iter(metric.keys())
        hull = plot.convex_hull(pos)

        for point in pos:
            min_distance = math.inf
            key = next(keys)
            for edge in hull:
                vector = plot.Vector.vec(point, edge)
                distance = plot.Vector.vec_len(vector)
                if distance < min_distance:
                    min_distance = distance
            if min_distance <= b_opt:
                # requires correction
                value = metric[key]
                delta = (m_opt * b_opt) - (m_opt * min_distance)
                corrected_metric[key] = value + delta
            else:
                # unaffected by boundary
                corrected_metric[key] = metric[key]

        return corrected_metric
