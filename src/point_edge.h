#pragma once

#include <iostream>
#include <fstream>


#include <boost/tuple/tuple.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// CGAL
#include <CGAL/number_utils.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   AK;
typedef CGAL::Polygon_2<AK>                           Polygon_2;

#include <utility> // defines std::pair

// #include "line_shape.cpp"
#include "region_growing.h"
// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Line_3 Line;
// Point with normal vector stored in a std::pair.
typedef boost::tuple<Point, Vector, int, bool, double, int, bool, double, int, bool> PNL;
typedef CGAL::Nth_of_tuple_property_map<0, PNL> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNL> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNL> Label_map;
typedef CGAL::Nth_of_tuple_property_map<3, PNL> IsWall_map;
typedef CGAL::Nth_of_tuple_property_map<4, PNL> LineFit_map;
typedef CGAL::Nth_of_tuple_property_map<5, PNL> JumpCount_map;
typedef CGAL::Nth_of_tuple_property_map<6, PNL> IsStep_map;
typedef CGAL::Nth_of_tuple_property_map<7, PNL> JumpEle_map;
typedef CGAL::Nth_of_tuple_property_map<8, PNL> Id_map;
typedef CGAL::Nth_of_tuple_property_map<9, PNL> IsHorizontal_map;
typedef std::vector<PNL>                        PNL_vector;
// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// search tree
// typedef boost::tuple<Point_3,int>                           Point_and_int;
// typedef CGAL::Random_points_in_cube_3<Point_3>              Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                       Traits_base;
typedef CGAL::Search_traits_adapter<PNL,
  CGAL::Nth_of_tuple_property_map<0, PNL>,
  Traits_base>                                              TreeTraits;
// typedef CGAL::Search_traits_3<SCK> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

// least squares stuff
// typedef CGAL::Simple_cartesian<double> SCK;
typedef CGAL::Cartesian<double> SCK;
typedef SCK::Point_3 Point_SCK;
typedef SCK::Line_3 Line_SCK;




namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> point_type;
typedef bg::model::point<double, 3, bg::cs::cartesian> point_type_3d;
typedef bg::model::segment<point_type> segment;

struct config {
  public:
  int metrics_normal_k = 10;
  int metrics_plane_min_points = 25;
  float metrics_plane_epsilon = 0.2;
  float metrics_plane_normal_threshold = 0.75;
  float metrics_is_wall_threshold = 0.3;
  float metrics_is_horizontal_threshold = 0.9;
  int metrics_k_linefit = 15;
  int metrics_k_jumpcnt_elediff = 10;

  int classify_jump_count_min = 1;
  int classify_jump_count_max = 5;
  float classify_line_dist = 0.005;
  float classify_jump_ele = 1.0;

  float linedetect_dist_threshold = 0.3;
  int linedetect_min_segment_count = 8;
  int linedetect_k = 10;

  float step_height_threshold = 1.0;
  float zrange_threshold = 0.2;
  bool merge_segid = true;
  bool merge_zrange = false;
  bool merge_step_height = true;
  bool merge_unsegmented = false;
  bool merge_dangling_egdes = false;
};


Polygon_2 ring_to_cgal_polygon(geoflow::LinearRing& ring);