%module ce_updater
%include "exception.i"
%include <std_string.i>
%include <std_map.i>
%include <std_vector.i>

%exception {
  try {
    $action
  } catch(const std::exception &e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#define SWIG_FILE_WITH_INIT

#include "ce_updater.hpp"
#include "wang_landau_sampler.hpp"
#include "adaptive_windows.hpp"
#include "cluster_tracker.hpp"
#include "pair_constraint.hpp"
%}
%include "numpy.i"

%init %{
  import_array();
%}

%template(map_str_dbl) std::map<std::string,double>;
%template(string_vector) std::vector<std::string>;

%include "ce_updater.hpp"
%include "matrix.hpp"
%include "matrix.tpp"
%include "cf_history_tracker.hpp"
%include "wang_landau_sampler.hpp"
%include "histogram.hpp"
%include "adaptive_windows.hpp"
%include "additional_tools.hpp"
%include "additional_tools.tpp"
%include "mc_observers.hpp"
%include "cluster_tracker.hpp"
%include "pair_constraint.hpp"
