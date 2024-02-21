#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <boost/heap/fibonacci_heap.hpp>

#include <geoflow/geoflow.hpp>

namespace tinsimp {

// fibonacci heap for greedy insertion code
struct point_error {
  point_error(size_t i, double e) : index(i), error(e){}
  size_t index;
  double error;
  size_t line_id;
  
  bool operator<(point_error const & rhs) const
  {
    return error < rhs.error;
  }
};
typedef boost::heap::fibonacci_heap<point_error> Heap;
typedef Heap::handle_type heap_handle;

typedef CGAL::Exact_predicates_inexact_constructions_kernel			K;
typedef CGAL::Projection_traits_xy_3<K>								Gt;
typedef CGAL::Triangulation_vertex_base_with_id_2<Gt>				Vb;
struct FaceInfo2
{
  FaceInfo2() {}
  std::vector<heap_handle>* points_inside = nullptr;
  CGAL::Plane_3<K>* plane = nullptr;
  int nesting_level = -1;
  bool in_domain() {
      return nesting_level % 2 == 1;
    }
};
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Gt>	Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Gt, Fbb>		Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>				Tds;
typedef CGAL::Exact_predicates_tag									Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>	CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>     CT;
typedef CDT::Point													Point;

// void greedy_insert(CDT &T, const std::vector<std::array<float,3>> &pts, double threshold);
void greedy_insert(CDT &T, std::vector<Point> &cpts, double threshold, size_t minpts=0);
void mark_domains(CDT& cdt);
}