#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Random.h>

#include <cstddef>
#include <gmp.h>
#include <vector>
#include <unordered_set>

#include "tinsimp.hpp"

namespace tinsimp {

void mark_domains(CDT& ct,
  CDT::Face_handle start,
  int index,
  std::list<CDT::Edge>& border) {
  if (start->info().nesting_level != -1) {
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while (!queue.empty()) {
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if (fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for (int i = 0; i < 3; i++) {
        CDT::Edge e(fh, i);
        CDT::Face_handle n = fh->neighbor(i);
        if (n->info().nesting_level == -1) {
          if (ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}
/**
* mark the triangles that are inside the original 2D polygon.
* explore set of facets connected with non constrained edges,
* and attribute to each such set a nesting level.
* start from facets incident to the infinite vertex, with a nesting
* level of 0. Then recursively consider the non-explored facets incident
* to constrained edges bounding the former set and increase the nesting level by 1.
* facets in the domain are those with an odd nesting level.
*/
void mark_domains(CDT& cdt) {
  // for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
  //   it->info().nesting_level = -1;
  // }
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while (!border.empty()) {
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if (n->info().nesting_level == -1) {
      mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
    }
  }
}

struct PointXYHash {
  std::size_t operator()(Point const& p) const noexcept {
    std::size_t h1 = std::hash<double>{}(p.x());
    std::size_t h2 = std::hash<double>{}(p.y());
    return h1 ^ (h2 << 1);
  }
};
struct PointXYEqual {
  bool operator()(Point const& p1, Point const& p2) const noexcept {
    auto ex = p1.x() == p2.x();
    auto ey = p1.y() == p2.y();
    return ex && ey;
  }
};

inline double compute_error(Point &p, CDT::Face_handle &face);

//--- TIN Simplification
// Greedy insertion/incremental refinement algorithm adapted from "Fast polygonal approximation of terrain and height fields" by Garland, Michael and Heckbert, Paul S.
inline double compute_error(Point &p, CDT::Face_handle &face) {
  if(!face->info().plane)
    face->info().plane = new CGAL::Plane_3<K>(
      face->vertex(0)->point(),
      face->vertex(1)->point(),
      face->vertex(2)->point());
  if(!face->info().points_inside)
    face->info().points_inside = new std::vector<heap_handle>();

  auto plane = face->info().plane;
  auto interpolate = -plane->a()/plane->c() * p.x() - plane->b()/plane->c()*p.y() - plane->d()/plane->c();
  double error = std::fabs(interpolate - p.z());
  return error;
}

// void greedy_insert(CDT &T, const std::vector<std::array<float,3>> &pts, double threshold) {
//   // Convert all elevation points to CGAL points
//   std::vector<Point> cpts;
//   cpts.reserve(pts.size());
//   for (auto& p : pts) {
//     cpts.push_back(Point(p[0], p[1], p[2]));
//   }

//   greedy_insert(T, cpts, threshold);
// }

void greedy_insert(CDT &T, std::vector<Point> &cpts, double threshold, size_t minpts) {
  // assumes all lidar points are inside a triangle
  Heap heap;

  // compute initial point errors, build heap, store point indices in triangles
  {
    std::unordered_set<Point, PointXYHash, PointXYEqual> set;
    for (int i = 0; i < cpts.size(); i++) {
      auto p = cpts[i];
      CDT::Locate_type lt;
      int li;
      auto not_duplicate = set.insert(p).second;
      if(not_duplicate){
        CDT::Face_handle face = T.locate(p, lt, li);
        if (lt == CDT::EDGE || lt == CDT::FACE) {
          auto e = compute_error(p, face);
          auto handle = heap.push(point_error(i, e));
          face->info().points_inside->push_back(handle);
        }
      }
    }
  }
  std::cout << "prepared tinsimp...\n";

  // insert points, update errors of affected triangles until threshold error is reached
  size_t insert_cnt = 0;
  // while ( (!heap.empty() && heap.top().error > threshold) || insert_cnt < minpts ) {
  while (!heap.empty() && heap.top().error > threshold) {
    // get top element (with largest error) from heap
    auto maxelement = heap.top();
    auto max_p = cpts[maxelement.index];

    // get triangles that will change after inserting this max_p
    std::vector<CDT::Face_handle> faces;
    T.get_conflicts(max_p, std::back_inserter(faces));

    // insert max_p in triangulation
    auto face_hint = faces[0];
    auto v = T.insert(max_p, face_hint);
    face_hint = v->face();
    ++insert_cnt;

    // update clear info of triangles that just changed, collect points that were inside these triangles
    std::vector<heap_handle> points_to_update;
    for (auto face : faces) {
      if (face->info().plane) {
        delete face->info().plane;
        face->info().plane = nullptr;
      }
      if (face->info().points_inside) {
        for (auto h : *face->info().points_inside) {
          if (maxelement.index != (*h).index)
            points_to_update.push_back(h);
        }
        std::vector<tinsimp::heap_handle>().swap((*face->info().points_inside));
      }
    }

    // remove the point we just inserted in the triangulation from the heap
    heap.pop();

    // update the errors of affected elevation points
    for (auto curelement : points_to_update) {
      auto p = cpts[(*curelement).index];
      //auto containing_face = T.locate(p, face_hint);
      CDT::Locate_type lt;
      int li;
      CDT::Face_handle containing_face = T.locate(p, lt, li, face_hint);
      if (lt == CDT::EDGE || lt == CDT::FACE) {
        const double e = compute_error(p, containing_face);
        const point_error new_pe = point_error((*curelement).index, e);
        heap.update(curelement, new_pe);
        containing_face->info().points_inside->push_back(curelement);
      }
    }
  }
  
  // insert more points so we can guarantee the minpts number of interior points
  if (heap.size() < minpts) minpts = heap.size();
  if(insert_cnt < minpts) {
    std::list<CGAL::Point_3<K>> remaining_pts;
    for (auto hit = heap.ordered_begin(); hit != heap.ordered_end(); ++hit) {
      auto& p = cpts[(*hit).index];
      remaining_pts.push_back(p);
    }

    auto rand = CGAL::Random(13374269);
    while (insert_cnt < minpts) {
      auto it = remaining_pts.begin();
      auto r = rand.get_int(0, remaining_pts.size());
      while (r > 0) {
        ++it;
        --r;
      }
      T.insert(*it);
      ++insert_cnt;
      remaining_pts.erase(it);
    }
  }

  //cleanup the stuff I put in face info of triangles
  for (CDT::Finite_faces_iterator fit = T.finite_faces_begin();
    fit != T.finite_faces_end(); ++fit) {
    if (fit->info().plane) {
      delete fit->info().plane;
      fit->info().plane = nullptr;
    }
    if (fit->info().points_inside) {
      delete fit->info().points_inside;
      fit->info().points_inside = nullptr;
    }
  }

}

}