// triangulation
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>

#include "stepedge_nodes.hpp"

namespace geoflow::nodes::stepedge {

  void Filter25DNode::process(){
    struct FaceInfo
    {
      FaceInfo(){}
      bool is_big=false;
    };
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Projection_traits_xy_3<K>                     Gt;
    typedef CGAL::Triangulation_vertex_base_2<Gt>               Vbb;
    typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>    Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Gt> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        TDS;
    typedef CGAL::Delaunay_triangulation_2<Gt, TDS>             DT_;
    typedef CGAL::Triangulation_hierarchy_2<DT_>                DT;
    typedef DT::Edge_circulator                                 Edge_circulator;

    DT dt;
    arr3f boxmin;
    arr3f boxmax;
    if(input("points").is_connected_type(typeid(PointCollection))) {
      auto& points = input("points").get<PointCollection&>();
      for (auto& p : points) {
        dt.insert(DT::Point(p[0], p[1], p[2]));
      }
      auto box = points.box();
      boxmin = box.min();
      boxmax = box.max();
    } else {
      auto& points_per_plane = input("points").get<IndexedPlanesWithPoints&>();
      Box box;
      for (auto& [plane_id, plane_pts] : points_per_plane) {
        if (plane_id<1) continue;
        for (auto& p : plane_pts.second) {
          dt.insert(DT::Point(p.x(), p.y(), p.z()));
          box.add({float(p.x()), float(p.y()), float(p.z())});
        }
      }
      boxmin = box.min();
      boxmax = box.max();
    }
    // build triangulation

    // mark low vertices
    size_t cnt_total, cnt_sharp;
    K::Vector_3 up(0.0,0.0,1.0);
    std::vector<DT::Vertex_handle> to_remove;
    PointCollection pts_removed, pts_remaining;
    for (auto& v: dt.finite_vertex_handles()) {
      cnt_total = cnt_sharp = 0;
      Edge_circulator ec = dt.incident_edges(v),
      done(ec);
      if (ec != 0) {
        do {
          K::Vector_3 a(v->point(), ec->first->vertex(dt.ccw(ec->second))->point());
          K::Vector_3 b(v->point(), ec->first->vertex(dt.cw(ec->second))->point());
          // auto angle = CGAL::approximate_angle(a,up);
          auto tri_angle = CGAL::approximate_angle(a,b);
          // auto up_angle = std::min(CGAL::approximate_angle(a,up), CGAL::approximate_angle(b,up));
          auto up_angle = CGAL::approximate_angle(a,up);
          // std::cout << angle << "\n";
          if ((up_angle < angle_thres)) {
            ++cnt_sharp;
          }
          ++cnt_total;
          // std::cout << dt.segment(ec) << std::endl;
          // compute polar angle of ec
        } while (++ec != done);
      }
      auto p = v->point();
      if (cnt_sharp >= count_thres) {
        to_remove.push_back(v);
        pts_removed.push_back({float(p.x()), float(p.y()), float(p.z())});
      } else {
        pts_remaining.push_back({float(p.x()), float(p.y()), float(p.z())});
      }
    }
    // delete low vertices
    for (auto& v : to_remove) {
      dt.remove(v);
    }

    // mark big triangles
    TriangleCollection triangles;
    auto sq_area_thres = area_thres*area_thres;
    for (auto& fh : dt.finite_face_handles()) {
      if( dt.triangle(fh).squared_area() > sq_area_thres ) {
        fh->info().is_big = true;
      } else {
        arr3f p0 = {float (fh->vertex(0)->point().x()), float (fh->vertex(0)->point().y()), float (fh->vertex(0)->point().z())};
        arr3f p1 = {float (fh->vertex(1)->point().x()), float (fh->vertex(1)->point().y()), float (fh->vertex(1)->point().z())};
        arr3f p2 = {float (fh->vertex(2)->point().x()), float (fh->vertex(2)->point().y()), float (fh->vertex(2)->point().z())
          };
        triangles.push_back({ p0,p1,p2 });
      }
    }

    // build heightfield
    RasterTools::Raster r(cellsize, boxmin[0], boxmax[0], boxmin[1], boxmax[1]);
    r.prefill_arrays(RasterTools::MAX);
    PointCollection heightfield_pts;
    for(size_t col=0; col<r.dimx_; ++col) {
      for(size_t row=0; row<r.dimy_; ++row) {
        auto p = r.getPointFromRasterCoords(col, row);
        // do linear TIN interpolation
        auto face = dt.locate(as::K::Point_3(double(p[0]),double(p[1]),0));
        if (face==nullptr) continue;
        if (dt.is_infinite(face)) continue;
        if (face->info().is_big) continue;
        CGAL::Plane_3<as::K> plane(
          face->vertex(0)->point(),
          face->vertex(1)->point(),
          face->vertex(2)->point()
        );
        double z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
        r.set_val(col, row, z_interpolate);
        p[2] = z_interpolate;
        heightfield_pts.push_back(p);
      }
    }

    output("points").set(pts_remaining);
    output("points_filtered").set(pts_removed);
    output("triangles").set(triangles);
    output("heightfield").set(r);
    output("heightfield_pts").set(heightfield_pts);
  }


  void PCRasteriseNode::process() {
    auto& points = input("points").get<PointCollection&>();
    auto box = points.box();
    auto boxmin = box.min();
    auto boxmax = box.max();

    RasterTools::Raster r(cellsize, boxmin[0], boxmax[0], boxmin[1], boxmax[1]);
    r.prefill_arrays(RasterTools::MAX);

    for(auto& p : points) {
      r.add_point(p[0], p[1], p[2], RasterTools::MAX);
    }
    PointCollection grid_points;
    vec1f values;
    double nodata = r.getNoDataVal();
    for(size_t i=0; i<r.dimx_ ; ++i) {
      for(size_t j=0; j<r.dimy_ ; ++j) {
        auto p = r.getPointFromRasterCoords(i,j);
        if (p[2]!=nodata) {
          grid_points.push_back(p);
          values.push_back(p[2]);
        }
      }
    }
    output("heightfield").set(r);
    output("values").set(values);
    output("grid_points").set(grid_points);
  }

  void SegmentRasteriseNode::process() {
    auto& alpha_rings = input("alpha_rings").get<LinearRingCollection&>();
    // auto& alpha_dts = input("alpha_dts").get<std::vector<as::Triangulation_2>&>();
    auto& roofplane_ids = input("roofplane_ids").get<vec1i&>();
    auto& pts_per_roofplane = input("pts_per_roofplane").get<IndexedPlanesWithPoints&>();
    RasterTools::Raster r;
    // if (input("heightfield").has_connection()) {
      // r = input("heightfield").get<RasterTools::Raster>();
    // } else {
      auto box = alpha_rings.box();
      auto boxmin = box.min();
      auto boxmax = box.max();
      r = RasterTools::Raster(cellsize, boxmin[0], boxmax[0], boxmin[1], boxmax[1]);
      r.prefill_arrays(RasterTools::MAX);
    // }

    PointCollection grid_points;
    size_t ring_cntr=0;
    for(auto& polygon : alpha_rings) {
      auto points_inside = r.rasterise_polygon(polygon);
      // auto dt = alpha_dts[ring_cntr];
      as::Triangulation_2 T;
      auto points = pts_per_roofplane[roofplane_ids[ring_cntr]].second;
      T.insert(points.begin(), points.end());
      auto as = as::Alpha_shape_2(T, as::FT(thres_alpha), as::Alpha_shape_2::GENERAL);

      for (auto& p : points_inside) {
        // grid_points.push_back(p);
        
        // do linear TIN interpolation
        auto face = as.locate(as::K::Point_3(double(p[0]),double(p[1]),0));
        if (face==nullptr) continue;
        if (as.is_infinite(face)) continue;
        if (as.classify(face) == as::Alpha_shape_2::EXTERIOR) continue;

        CGAL::Plane_3<as::K> plane(
          face->vertex(0)->point(),
          face->vertex(1)->point(),
          face->vertex(2)->point()
        );
        double z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
        r.add_point(p[0], p[1], z_interpolate, RasterTools::MAX);

        // do plane projection
        // auto& plane = pts_per_roofplane[roofplane_ids[ring_cntr]].first;
        // double z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
        // r.add_point(p[0], p[1], z_interpolate, RasterTools::MAX);
      }
      ++ring_cntr;
    }
    vec1f values;
    double nodata = r.getNoDataVal();
    for(size_t i=0; i<r.dimx_ ; ++i) {
      for(size_t j=0; j<r.dimy_ ; ++j) {
        auto p = r.getPointFromRasterCoords(i,j);
        if (p[2]!=nodata) {
          grid_points.push_back(p);
          values.push_back(p[2]);
        }
      }
    }
    output("heightfield").set(r);
    output("values").set(values);
    output("grid_points").set(grid_points);
  }

  // void GridMaxNode::process() {
  //   auto& g1 = input("grid_1").get<RasterTools::Raster>();
  //   auto& g2 = input("grid_2").get<RasterTools::Raster>();

  //   RasterTools::Raster r(g1.cellSize_, g1.minx_, g1.maxx_, g1.miny_, g1.maxy_);
  //   r.prefill_arrays(RasterTools::MAX);

  //   PointCollection grid_points;
  //   for(size_t i=0; i<r.dimx_ ; ++i) {
  //     for(size_t j=0; j<r.dimy_ ; ++j) {
  //       size_t l = j*r.dimx_+i;
  //       r.vals_[l] = std::max(g1.vals_[l], g2.vals_[l]);
  //       auto p = r.getPointFromRasterCoords(i,j);
  //       if(p[2]!= r.noDataVal_)
  //         grid_points.push_back(p);
  //     }
  //   }

  //   output("heightfield").set(r);
  //   output("grid_points").set(grid_points);
  // }
}