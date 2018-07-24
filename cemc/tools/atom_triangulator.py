from scipy.spatial import Delaunay
import numpy as np


class AtomTriangulator(object):
    """Class that creates a Delaunay triangulation of a set of atoms."""
    def __init__(self, atoms):
        self.atoms = atoms

    def triangulate(self):
        """Performs a Delaunay triangulation."""
        return Delaunay(self.atoms.get_positions())

    def get_volume(self, points, simplex):
        """Compute the volume of one simplex.

        :param points: Nx3 array of the points of the atoms
        :param simplex: One simplex from the Delaunay object
        """
        vectors = np.zeros((3, 3))
        vectors[:, 0] = points[simplex[1], :] - points[simplex[0], :]
        vectors[:, 1] = points[simplex[2], :] - points[simplex[0], :]
        vectors[:, 2] = points[simplex[3], :] - points[simplex[0], :]
        return np.det(vectors) / 6.0

    def get_volumes(self, delaunay):
        """Compute the volume of all tetrahedons.

        :param delaunay: Instance of the Delaunay object
        """
        points = self.atoms.get_positions()

        volumes = []
        for simplex in delanay.simplices:
            V = self.get_volume(points, simplex)
            volumes.append(V)
        return volumes

    def get_volume_statistics(self, delaunay):
        """Compute statistics of the volume of the trianges.

        :param delaunay: Instance of the Delaunay object
        """
        points = self.atoms.get_positions()

        volumes = self.get_volumes(delaunay)
        for simplex in delaunay.simplices:
            V = self.get_volume(points, simplex)
            volumes.append(V)

        statistics = {
            "max": np.max(volumes),
            "min": np.min(volumes),
            "average": np.mean(volumes),
            "std": np.std(volumes),
            "hist": np.histogram(volumes, bins="auto")
        }
        return statistics

    def filter_volume_percentile(self, delaunay, percentile=None,
                                 max_volume=None):
        """Filter simplex list based on percentile.

        :param delaunay: Instance of the Delaunay object
        :param percentile: If given, keep the smallest.
            If percentils is 95, the 95 percent smallest clusters will be kept.
            Note that if both percentile and max_volume is given, max_volume
            will be prioritized.
        :param max_volume: If given Maximum volume to keep
        """
        volumes = self.get_volumes(delaunay)

        if max_volume is not None:
            pass
        elif percentile is not None:
            max_volume = np.percentile(volumes, percentile)
        else:
            raise ValueError("Have to give either max_volume or percentile!")

        filtered_simplices = []
        points = self.atoms.get_positions()
        for simplex in delaunay.simplices:
            V = self.get_volume(points, simplex)
            if V < max_volume:
                filtered_simplices.append(simplex)
        delaunay.simplices = filtered_simplices
        return delaunay
