# distutils: language = c++

from cemc.cpp_ext.ce_updater cimport CEUpdater

cdef class PyCEUpdater:
    """
    Cython wrapper for the C++ class
    """
    cdef CEUpdater _cpp_class

    def __init__(self, bc, corr_func, eci):
        self._cpp_class.init(bc, corr_func, eci)

    def clear_history(self):
        self._cpp_class.clear_history()

    def undo_changes(self):
        self._cpp_class.undo_changes()

    def update_cf(self, system_changes):
        self._cpp_class.update_cf(system_changes)

    def add_linear_vib_correction(self, value):
        self._cpp_class.add_linear_vib_correction(value)

    def vid_energy(self, T):
        return self._cpp_class.vib_energy(T)

    def get_cf(self):
        return self._cpp_class.get_cf()

    def set_ecis(self, ecis):
        self._cpp_class.set_ecis(ecis)

    def get_singlets(self):
        return self._cpp_class.get_singlets()

    def get_singlets(self, array):
        self._cpp_class.get_singlets(array)
