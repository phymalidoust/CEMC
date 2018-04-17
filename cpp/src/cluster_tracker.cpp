#include "cluster_tracker.hpp"
#include "cluster.hpp"
#include "matrix.hpp"
#include "additional_tools.hpp"
#include <stdexcept>
#include <sstream>

using namespace std;

ClusterTracker::ClusterTracker( CEUpdater &updater, const std::string &cname, const std::string &element ): \
updater(&updater),cname(cname),element(element)
{
  verify_cluster_name_exists();
};

void ClusterTracker::find_clusters()
{
  const vector<string> &symbs = updater->get_symbols();
  atomic_clusters.resize( symbs.size() );
  for ( unsigned int i=0;i<atomic_clusters.size();i++ )
  {
    atomic_clusters[i] = -1; // All atoms are initially a root site
  }
  const vector< map<string,Cluster> >& clusters = updater->get_clusters();
  const Matrix<int>& trans_mat = updater->get_trans_matrix();

  for ( unsigned int i=0;i<symbs.size();i++ )
  {
    // If the element does not match, do not do anything
    if ( symbs[i] != element )
    {
      continue;
    }

    int current_root_indx = i;
    while( atomic_clusters[current_root_indx] != -1 )
    {
      current_root_indx = atomic_clusters[current_root_indx];
    }

    if ( atomic_clusters[current_root_indx] != -1 )
    {
      throw runtime_error( "Something strange happend. ID of root index is not -1!" );
    }

    // Loop over all symmetries
    for ( unsigned int trans_group=0;trans_group<clusters.size();trans_group++ )
    {
      if ( clusters[trans_group].find(cname) == clusters[trans_group].end() )
      {
        // Cluster does not exist in this translattional symmetry group
        continue;
      }

      const std::vector< std::vector<int> >& members = clusters[trans_group].at(cname).get();
      for ( int subgroup=0;subgroup<members.size();subgroup++ )
      {
        int indx = trans_mat(i,members[subgroup][0]);

        if ( symbs[indx] == element )
        {
          int root_indx = indx;
          while( atomic_clusters[root_indx] != -1 )
          {
            root_indx = atomic_clusters[root_indx];
          }

          if ( root_indx != current_root_indx )
          {
            atomic_clusters[root_indx] = current_root_indx;
          }
        }
      }
    }
  }
}

void ClusterTracker::get_cluster_statistics( map<string,double> &res, vector<int> &cluster_sizes ) const
{
  double average_size = 0.0;
  double max_size = 0.0;
  double avg_size_sq = 0.0;
  map<int,int> num_members_in_cluster;
  cluster_sizes.clear();

  for ( unsigned int i=0;i<atomic_clusters.size();i++ )
  {
    int root_indx = i;
    while ( atomic_clusters[root_indx] != -1 )
    {
      root_indx = atomic_clusters[root_indx];
    }

    if ( root_indx != i )
    {
      if ( num_members_in_cluster.find(root_indx) != num_members_in_cluster.end() )
      {
        num_members_in_cluster[root_indx] += 1;
      }
      else
      {
        num_members_in_cluster[root_indx] = 1;
      }
    }
  }

  for ( auto iter=num_members_in_cluster.begin(); iter != num_members_in_cluster.end(); ++iter )
  {
    if ( iter->second > 2 )
    {
      cluster_sizes.push_back(iter->second);
    }
    average_size += iter->second;
    avg_size_sq += iter->second*iter->second;
    if ( iter->second > max_size )
    {
      max_size = iter->second;
    }
  }
  res["avg_size"] = average_size;
  res["max_size"] = max_size;
  res["avg_size_sq"] = avg_size_sq;
  res["number_of_clusters"] = cluster_sizes.size();
}

PyObject* ClusterTracker::get_cluster_statistics_python() const
{
  PyObject* dict = PyDict_New();
  map<string,double> res;
  vector<int> cluster_sizes;
  get_cluster_statistics(res,cluster_sizes);
  for ( auto iter=res.begin(); iter != res.end(); ++iter )
  {
      PyObject* value = PyFloat_FromDouble( iter->second );
      PyDict_SetItemString( dict, iter->first.c_str(), value );
      Py_DECREF(value);
  }

  PyObject* size_list = PyList_New(0);
  for ( int i=0; i< cluster_sizes.size();i++ )
  {
    PyObject *value = PyInt_FromLong( cluster_sizes[i] );
    PyList_Append( size_list, value );
    Py_DECREF(value);
  }
  PyDict_SetItemString( dict, "cluster_sizes", size_list );

  return dict;
}

void ClusterTracker::atomic_clusters2group_indx( vector<int> &group_indx ) const
{
  group_indx.resize( atomic_clusters.size() );
  for ( unsigned i=0;i<atomic_clusters.size();i++ )
  {
    int root_indx = i;
    while ( atomic_clusters[root_indx] != -1 )
    {
      root_indx = atomic_clusters[root_indx];
    }
    group_indx[i] = root_indx;
  }
}

PyObject* ClusterTracker::atomic_clusters2group_indx_python() const
{
  PyObject *list = PyList_New(0);
  vector<int> grp_indx;
  atomic_clusters2group_indx(grp_indx);
  for ( unsigned int i=0;i<grp_indx.size();i++ )
  {
    PyObject *pyint = PyInt_FromLong( grp_indx[i] );
    PyList_Append( list, pyint );
    Py_DECREF(pyint);
  }
  return list;
}

void ClusterTracker::verify_cluster_name_exists() const
{
  const vector< map<string,Cluster> >& clusters = updater->get_clusters();
  vector<string> all_names;
  for ( unsigned int i=0;i<clusters.size();i++ )
  {
    if ( clusters[i].find(cname) != clusters[i].end() )
    {
      return;
    }
    for ( auto iter=clusters[0].begin(); iter != clusters[0].end(); ++iter )
    {
      all_names.push_back( iter->first );
    }
  }
  stringstream ss;
  ss << "There are now correlation functions corresponding to the cluster name given!\n";
  ss << "Available names:\n";
  ss << all_names;
  throw invalid_argument( ss.str() );
}
