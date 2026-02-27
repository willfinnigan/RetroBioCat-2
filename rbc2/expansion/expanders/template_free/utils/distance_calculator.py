"""Pure Python/scipy replacement for the C extension BFS distance calculator.

The original ASKCOS code uses a compiled C extension (c_calculate.cpp/.so) for
BFS all-pairs shortest path on molecular graphs. For single-molecule inference
(n=20-80 nodes), scipy's BFS shortest path is trivially fast.
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path


class DistanceCalculator:

    @staticmethod
    def calculate(adjacency: np.ndarray,
                  a_length: int,
                  max_distance: int) -> np.ndarray:
        """Calculate all-pairs shortest path distances on a graph.

        Args:
            adjacency: Boolean adjacency matrix of shape (a_length, a_length).
            a_length: Number of nodes.
            max_distance: Maximum distance to consider.

        Returns:
            Distance matrix of shape (a_length, a_length) with dtype uint8.
        """
        sparse_adj = csr_matrix(adjacency.astype(np.float32))
        dist_matrix = shortest_path(sparse_adj, method='D', unweighted=True)

        # Unreachable nodes (inf) get max_distance, then clamp and convert
        dist_matrix[np.isinf(dist_matrix)] = max_distance
        dist_matrix = np.clip(dist_matrix, 0, max_distance).astype(np.uint8)

        return dist_matrix
