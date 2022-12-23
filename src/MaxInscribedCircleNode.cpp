#include "MaxInscribedCircleNode.hpp"

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
// #include <CGAL/Voronoi_diagram_2.h>
// #include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
// #include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

#include <CGAL/squared_distance_2.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/oriented_side.h>
#include <CGAL/Polygon_2.h> 
// #include <CGAL/Polygon_2_algorithms.h>

#include <chrono>

#include "ptinpoly.h"

static const double PI = 3.141592653589793238462643383279502884;

namespace geoflow::nodes::stepedge {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Delaunay_triangulation_2<K>   Triangulation;
  typedef Triangulation::Edge_iterator        Edge_iterator;
  typedef Triangulation::Point                Point;
  typedef CGAL::Polygon_2<K>                  Polygon;
  typedef CGAL::Polygon_with_holes_2<K>       Polygon_with_holes;

  void insert_edges(Triangulation& t, const Polygon& polygon, const float& interval) {
    for (auto ei = polygon.edges_begin(); ei != polygon.edges_end(); ++ei) {
      auto e_l = CGAL::sqrt(ei->squared_length());
      auto e_v = ei->to_vector();
      auto n = std::floor(e_l / interval);
      auto s = ei->source();
      t.insert(s);
      for(size_t i=0; i<n; ++i) {
        t.insert(s + i*interval*e_v);
        // l += interval;
      }
    }
  }

  class GridPIPTester {

    pGridSet ext_gridset;
    std::vector<pGridSet> hole_gridsets;
    int Grid_Resolution = 20;

    pGridSet build_grid(const Polygon& ring) {

      int size = ring.size();
      std::vector<pPipoint> pgon;
      for (auto pi = ring.vertices_begin(); pi != ring.vertices_end(); ++pi) {
        pgon.push_back(new Pipoint{ pi->x(),pi->y() });
      }
      pGridSet grid_set = new GridSet();
      // skip last point in the ring, ie the repetition of the first vertex
      GridSetup(&pgon[0], pgon.size(), Grid_Resolution, grid_set);
      for (int i = 0; i < size; i++) {
        delete pgon[i];
      }
      return grid_set;
    }

    public:
    GridPIPTester(const Polygon_with_holes& polygon) {
      ext_gridset = build_grid(polygon.outer_boundary());
      for (auto& hole : polygon.holes()) {
        hole_gridsets.push_back(build_grid(hole));
      }
    }
    ~GridPIPTester() {
      delete ext_gridset;
      for (auto& h : hole_gridsets) {
        delete h;
      }
    }

    bool test(const Point& p) {
      pPipoint pipoint = new Pipoint{p.x(),p.y()};
      bool inside = GridTest(ext_gridset, pipoint);
      if (inside) {
        for (auto& hole_gridset : hole_gridsets) {
          inside = inside && !GridTest(hole_gridset, pipoint);
          if (!inside) break;
        }
      }
      delete pipoint;
      return inside;
    }

  };

  
  void MaxInscribedCircleNode::process() {   
    std::clock_t c_start = std::clock(); // CPU time
 
    // build grid

    // build VD/DT
    Triangulation t;
    
    // insert points point cloud
    auto& pointcloud = input("pointcloud").get<PointCollection&>();
    for (auto& p : pointcloud) {
      t.insert(Point(p[0], p[1]));
    }
    
    // insert pts on footprint boundary
    auto& lr = input("polygon").get<LinearRing>();

    Polygon poly2;
    for (auto& p : lr) {
      poly2.push_back(Point(p[0], p[1]));
    }
    std::vector<Polygon> holes;
    for (auto& lr_hole : lr.interior_rings()) {
      Polygon hole;
      for (auto& p : lr_hole) {
        hole.push_back(Point(p[0], p[1]));
      }
      holes.push_back(hole);
    }
    auto polygon = Polygon_with_holes(poly2, holes.begin(), holes.end());

    // double l = 0;
    insert_edges(t, polygon.outer_boundary(), polygon_densify_);
    for (auto& hole : polygon.holes()) {
      insert_edges(t, hole, polygon_densify_);
    }
    // build gridset for point in polygon checks
    auto pip_tester = GridPIPTester(polygon);

    // std::cout << 1000.0 * (std::clock()-c_start) / CLOCKS_PER_SEC << "ms 1\n";

    // find VD node with largest circle
    double r_max = 0;
    Point c_max;
    for(const auto& face : t.finite_face_handles()) {
      // get the voronoi node
      auto c = t.dual(face);
      // check it is inside footprint polygon
      if(pip_tester.test(c)) {
        for(size_t i=0; i<3; ++i) {
          auto r = CGAL::squared_distance(c, face->vertex(i)->point());
          if (r>r_max) {
            r_max = r;
            c_max = c;
          }
        }
      }
    }
    r_max = std::sqrt(r_max);
    // std::cout << 1000.0 * (std::clock()-c_start) / CLOCKS_PER_SEC << "ms 1\n";

    // std::cout << "Max radius: " << r_max << std::endl;
    // std::cout << "Max radius center: " << c_max << std::endl;

    LinearRing circle;
    const double angle_step = PI/5;
    for (float a=0; a<2*PI; a+=angle_step) {
      circle.push_back({
        float(c_max.x() + std::cos(a)*r_max), 
        float(c_max.y() + std::sin(a)*r_max),
        0
      });
    }
    // PointCollection mc; mc.push_back({float(c_max.x()), float(c_max.y()), 0});
    output("max_circle").set(circle);
    output("max_diameter").set(float(2*r_max));
  }

}