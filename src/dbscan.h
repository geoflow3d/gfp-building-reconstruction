// This file is part of gfp-building-reconstruction
// Copyright (C) 2018-2022 Ravi Peters

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include <unordered_map>

#include <CGAL/property_map.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>


namespace dbscan {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel cgal_kernel;
  typedef cgal_kernel::Point_3 Point;
  typedef cgal_kernel::Line_3 Line;

  using namespace std;

  class dbscanClusterer {
    
    typedef std::pair<Point,size_t> point_index;
    typedef CGAL::Search_traits_3<cgal_kernel>                       Traits_base;
    typedef CGAL::Search_traits_adapter<point_index,
    CGAL::First_of_pair_property_map<point_index>,
    Traits_base>                                              TreeTraits;
    typedef CGAL::Kd_tree<TreeTraits> Tree;
    typedef CGAL::Fuzzy_sphere<TreeTraits> Fuzzy_sphere;

    vector<point_index> indexed_points;
    Tree tree;
    vector<bool> point_seed_flags;
    size_t region_counter=1;
    
    public:
    vector<size_t> point_segment_idx; // 0=unsegmented, maybe put this on the heap...
    unordered_map<size_t, Line> segment_shapes;
    double dist_thres = 0.2*0.2; //epsilon
    size_t min_segment_count = 20;

    LineDetector(vector<Point> &points);
    vector<size_t> get_point_indices(size_t shape_id);
    void detect();

    private:
    void grow_region(size_t seed_idx);
  };
}