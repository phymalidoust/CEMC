from scipy.spatial import Delaunay
import numpy as np


class AtomTriangulator(object):
    """Class that creates a Delaunay triangulation of a set of atoms."""
    def __init__(self, atoms):
        self.atoms = atoms
        self._triangulation = self._triangulate()

    def _triangulate(self):
        """Performs a Delaunay triangulation."""
        return Delaunay(self.atoms.get_positions())

    def _log(self, msg):
        """Prints log messages."""
        print(msg)

    def get_volume(self, points, simplex):
        """Compute the volume of one simplex.

        :param points: Nx3 array of the points of the atoms
        :param simplex: One simplex from the Delaunay object
        """
        vectors = np.zeros((3, 3))
        vectors[:, 0] = points[simplex[1], :] - points[simplex[0], :]
        vectors[:, 1] = points[simplex[2], :] - points[simplex[0], :]
        vectors[:, 2] = points[simplex[3], :] - points[simplex[0], :]
        return np.abs(np.linalg.det(vectors)) / 6.0

    def get_volumes(self):
        """Compute the volume of all tetrahedons.
        """
        points = self.atoms.get_positions()

        volumes = []
        for simplex in self._triangulation.simplices:
            V = self.get_volume(points, simplex)
            volumes.append(V)
        return volumes

    def get_volume_statistics(self):
        """Compute statistics of the volume of the trianges.
        """
        points = self.atoms.get_positions()

        volumes = self.get_volumes()
        for simplex in self._triangulation.simplices:
            V = self.get_volume(points, simplex)
            volumes.append(V)

        lengths = np.array(self.get_edge_lengths())

        statistics = {
            "max_vol": np.max(volumes),
            "min_vol": np.min(volumes),
            "average_vol": np.mean(volumes),
            "std_vol": np.std(volumes),
            "max_length": np.max(lengths),
            "min_length": np.min(lengths),
            "average_length": np.mean(lengths),
            "median_length": np.median(lengths),
            "std_length": np.std(lengths)
            # "hist": np.histogram(volumes, bins="auto")
        }
        return statistics

    def filter_volume(self, percentile=None,
                      max_volume=None):
        """Filter simplex list based on percentile.

        :param percentile: If given, keep the smallest.
            If percentils is 95, the 95 percent smallest clusters will be kept.
            Note that if both percentile and max_volume is given, max_volume
            will be prioritized.
        :param max_volume: If given Maximum volume to keep
        """
        volumes = self.get_volumes()

        if max_volume is not None:
            pass
        elif percentile is not None:
            max_volume = np.percentile(volumes, percentile)
        else:
            raise ValueError("Have to give either max_volume or percentile!")

        filtered_simplices = []
        points = self.atoms.get_positions()
        n_init = len(self._triangulation.simplices)
        for simplex in self._triangulation.simplices:
            V = self.get_volume(points, simplex)
            if V < max_volume:
                filtered_simplices.append(simplex)
        self._triangulation.simplices = filtered_simplices
        n_final = len(filtered_simplices)
        self._log("Removed {} simplices".format(n_init-n_final))

    def export_to_gmsh(self, fname):
        """Export the mesh to GMSH."""
        with open(fname, 'w') as out:
            # Export mesh format
            out.write("$MeshFormat\n")
            out.write("2.2 0 8\n")
            out.write("$EndMeshFormat\n")
            # Export node positions
            out.write("$Nodes\n")
            out.write("{}\n".format(len(self.atoms)))
            for atom in self.atoms:
                pos = atom.position
                out.write("{} {} {} {}\n".format(
                    atom.index, pos[0], pos[1], pos[2]))
            out.write("$EndNodes\n")

            # Export elements
            n_elem = len(self._triangulation.simplices)
            out.write("$Elements\n")
            out.write("{}\n".format(n_elem))
            elem_type = 4  # Four node tetre hedron
            n_tags = 0
            for elm_num, s in enumerate(self._triangulation.simplices):
                out.write("{} {} {} {} {} {} {}\n".format(
                    elm_num, elem_type, n_tags, s[0], s[1], s[2], s[3]))
            out.write("$EndElements\n")
        self._log("Mesh file written to {}".format(fname))

    def get_edge_lengths(self):
        """Compute the lengths of all the edges."""
        pos = self.atoms.get_positions()
        all_lengths = []
        for s in self._triangulation.simplices:
            lengths = []
            for i in range(0, 4):
                for j in range(i+1, 4):
                    vec = pos[s[j], :] - pos[s[i], :]
                    lengths.append(np.sqrt(vec.dot(vec)))
            all_lengths.append(lengths)
        return all_lengths

    def volume_to_edge_length_ratio(self):
        """Compute the volume to edge length ratio."""
        vols = np.array(self.get_volumes())
        lengths = self.get_edge_lengths()
        sum_lengths = np.zeros_like(vols)
        for i, l in enumerate(lengths):
            sum_lengths[i] = np.sum(np.array(l)**2)
        return vols**(2.0/3.0) / sum_lengths

    def max_to_min_length_ratio(self):
        """Compute the max to min ratio of the edge lengths."""
        lengths = self.get_edge_lengths()
        ratios = []
        for l in lengths:
            ratios.append(np.min(l)/np.max(l))
        return np.array(ratios)

    def filter_min_over_max_length(self, min_ratio=0.25):
        """Filter simplices based on the ratio between min and max length.abs

        Keep only simplices where min/max >= min_ratio
        """
        ratios = self.max_to_min_length_ratio()
        filtered = []
        init_num = len(self._triangulation.simplices)
        for ratio, s in zip(ratios, self._triangulation.simplices):
            if ratio > min_ratio:
                filtered.append(s)
        final_num = len(filtered)
        self._triangulation.simplices = filtered
        self._log("Ratio filter removed {} simplices".format(
            init_num - final_num))

    def filter_max_length(self, max_length=None):
        """Filter the simplices based on the maximum length.

        Keep all simplices smaller than max_length
        :param max_length: Maximum length fo keep
        """
        lengths = self.get_edge_lengths()
        filtered = []
        n_initial = len(self._triangulation.simplices)
        for l, s in zip(lengths, self._triangulation.simplices):
            if np.max(l) < max_length:
                filtered.append(s)
        n_final = len(filtered)
        self._triangulation.simplices = filtered
        self._log("Max length removed: {}".format(n_initial-n_final))
