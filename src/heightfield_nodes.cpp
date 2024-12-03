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
// triangulation
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Cartesian.h>
#include <cstddef>
#include <cmath>
#include <geoflow/common.hpp>

#include "Raster.h"
#include "point_edge.h"
#include "stepedge_nodes.hpp"
#include "pip_util.hpp"

namespace geoflow::nodes::stepedge {

  class PointCloud25DTriangulator {
    
    public: 
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    private:
    struct FaceInfo
    {
      FaceInfo(){}
      bool is_big=false;
    };
    typedef CGAL::Projection_traits_xy_3<K>                     Gt;
    typedef CGAL::Triangulation_vertex_base_2<Gt>               Vbb;
    typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>    Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Gt> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        TDS;
    typedef CGAL::Delaunay_triangulation_2<Gt, TDS>             DT_;
    typedef CGAL::Triangulation_hierarchy_2<DT_>                DT;
    typedef DT::Edge_circulator                                 Edge_circulator;
    typedef DT::Face_circulator                                 Face_circulator;

    DT dt;
    bool do_area_thres;
    float area_thres;
    bool do_angle_thres;
    float angle_thres;
    bool do_normal_angle_thres;
    float normal_angle_thres;
    bool do_len_thres;
    float len_thres;

    public:

    PointCloud25DTriangulator(
      bool do_area_thres,
      float area_thres,
      bool do_angle_thres,
      float angle_thres,
      bool do_normal_angle_thres,
      float normal_angle_thres,
      bool do_len_thres,
      float len_thres
    ) : do_area_thres(do_area_thres),
        area_thres(area_thres),
        do_angle_thres(do_angle_thres),
        angle_thres(angle_thres),
        do_normal_angle_thres(do_normal_angle_thres),
        normal_angle_thres(normal_angle_thres), 
        do_len_thres(do_len_thres),
        len_thres(len_thres)
    {};

    void insert_point(std::array<float,3>& p) {
      dt.insert(DT::Point(p[0], p[1], p[2]));
    }

    void insert_point(float x, float y, float z) {
      dt.insert(DT::Point(x,y,z));
    }

    bool get_z(std::array<float,3>& p, float& z_interpolate) {
      auto face = dt.locate(K::Point_3(double(p[0]),double(p[1]),0));
      if (face==nullptr) return false;
      if (dt.is_infinite(face)) return false;
      if (face->info().is_big) return false;
      CGAL::Plane_3<K> plane(
        face->vertex(0)->point(),
        face->vertex(1)->point(),
        face->vertex(2)->point()
      );
      z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
      return true;
    }

    void mark_big_triangles() {
      typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        // mark big triangles
      auto sq_area_thres = area_thres*area_thres;
      auto sq_len_thres = len_thres*len_thres;
      K::Vector_3 up(0.0,0.0,1.0);
      for (auto& fh : dt.finite_face_handles()) {
        fh->info().is_big = false;
        if (do_area_thres) {
          if( dt.triangle(fh).squared_area() > sq_area_thres ) {
            fh->info().is_big = fh->info().is_big || true;
          }
        }
        if (do_angle_thres) {
          if( 
            CGAL::approximate_angle(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point()) < angle_thres ||
            CGAL::approximate_angle(fh->vertex(1)->point(), fh->vertex(2)->point(), fh->vertex(0)->point()) < angle_thres ||
            CGAL::approximate_angle(fh->vertex(2)->point(), fh->vertex(0)->point(), fh->vertex(1)->point()) < angle_thres
          ) {
            fh->info().is_big = fh->info().is_big || true;
          }
        }
        if (do_normal_angle_thres) {
          if( 
            CGAL::approximate_angle(dt.triangle(fh).supporting_plane().orthogonal_vector(), up) > normal_angle_thres
          ) {
            fh->info().is_big = fh->info().is_big || true;
          }
        }
        if (do_len_thres)  {
          if( (K::Segment_3(fh->vertex(0)->point(), fh->vertex(1)->point()).squared_length() > sq_len_thres) ||  
              (K::Segment_3(fh->vertex(1)->point(), fh->vertex(2)->point()).squared_length() > sq_len_thres) || 
              (K::Segment_3(fh->vertex(2)->point(), fh->vertex(0)->point()).squared_length() > sq_len_thres) ) {
            fh->info().is_big = fh->info().is_big || true;
          }
        }
      }

    }

    void remove_marked() {
      std::vector<DT::Vertex_handle> to_remove;
      // PointCollection pts_removed, pts_remaining;
      for (auto& v: dt.finite_vertex_handles()) {
        bool remove_vertex = true;
        Face_circulator fc = dt.incident_faces(v),
        done(fc);
        if (fc != 0) {
          do {
            remove_vertex = remove_vertex && fc->info().is_big;
          } while (++fc != done);
        }
        auto p = v->point();
        if (remove_vertex) {
          to_remove.push_back(v);
          // pts_removed.push_back({float(p.x()), float(p.y()), float(p.z())});
        } else {
          // pts_remaining.push_back({float(p.x()), float(p.y()), float(p.z())});
        }
      }
      // delete outlier vertices
      for (auto& v : to_remove) {
        dt.remove(v);
      }
    }

    void get_triangles(TriangleCollection& triangles) {
      for (auto& fh : dt.finite_face_handles()) {
        if ( !fh->info().is_big ) {
          arr3f p0 = {float (fh->vertex(0)->point().x()), float (fh->vertex(0)->point().y()), float (fh->vertex(0)->point().z())};
          arr3f p1 = {float (fh->vertex(1)->point().x()), float (fh->vertex(1)->point().y()), float (fh->vertex(1)->point().z())};
          arr3f p2 = {float (fh->vertex(2)->point().x()), float (fh->vertex(2)->point().y()), float (fh->vertex(2)->point().z())
            };
          triangles.push_back({ p0,p1,p2 });
        }
      }
    }

  };

  void Filter25DNode::process(){
    
    // build triangulation
    PointCloud25DTriangulator pdt(
      do_area_thres,
      area_thres,
      do_angle_thres,
      angle_thres,
      do_normal_angle_thres,
      normal_angle_thres,
      do_len_thres,
      len_thres
    );
    arr3f boxmin;
    arr3f boxmax;
    if(input("points").is_connected_type(typeid(PointCollection))) {
      auto& points = input("points").get<PointCollection&>();
      for (auto& p : points) {
        pdt.insert_point(p);
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
          pdt.insert_point(p.x(), p.y(), p.z());
          box.add(arr3f{float(p.x()), float(p.y()), float(p.z())});
        }
      }
      boxmin = box.min();
      boxmax = box.max();
    }

    pdt.mark_big_triangles();

    // remove vertices with all marked neighbour triangles
    
    pdt.remove_marked();

    pdt.mark_big_triangles();

    // get triangles
    TriangleCollection triangles;
    pdt.get_triangles(triangles);

    // build heightfield
    RasterTools::Raster r(cellsize, boxmin[0], boxmax[0], boxmin[1], boxmax[1]);
    r.prefill_arrays(RasterTools::MAX);
    PointCollection heightfield_pts;
    float z_interpolate;
    for(size_t col=0; col<r.dimx_; ++col) {
      for(size_t row=0; row<r.dimy_; ++row) {
        auto p = r.getPointFromRasterCoords(col, row);
        // do linear TIN interpolation
        if( pdt.get_z(p, z_interpolate) ) {
          r.set_val(col, row, z_interpolate);
          p[2] = z_interpolate;
          heightfield_pts.push_back(p);
        }
      }
    }

    // output("points").set(pts_remaining);
    // output("points_filtered").set(pts_removed);
    output("triangles").set(triangles);
    output("heightfield").set(r);
    output("heightfield_pts").set(heightfield_pts);
  }


  void PCRasteriseNode::process() {
    auto& points = input("points").get<PointCollection&>();
    auto& footprint = input("footprint").get<LinearRing&>();

    Box box;
    if (use_footprint_) {
      box = footprint.box();
    } else {
      box = points.box();
    }
    auto boxmin = box.min();
    auto boxmax = box.max();

    RasterTools::Raster r_max(cellsize, boxmin[0], boxmax[0], boxmin[1], boxmax[1]);
    r_max.prefill_arrays(RasterTools::MAX);

    RasterTools::Raster r_min(r_max), r_fp(r_max);
    r_min.prefill_arrays(RasterTools::MIN);
    r_fp.prefill_arrays(RasterTools::MAX);

    std::vector<std::vector<float>> buckets(r_max.dimx_*r_max.dimy_);

    if (use_footprint_) {
      auto exterior = build_grid(footprint);
      std::vector<pGridSet> holes;
      for (auto& hole : footprint.interior_rings()) {
        holes.push_back(build_grid(hole));
      }
      
      for (size_t col = 0; col < r_fp.dimx_; ++col) {
        for (size_t row = 0; row < r_fp.dimy_; ++row) {
          auto p = r_fp.getPointFromRasterCoords(col, row);
          pPipoint pipoint = new Pipoint{p[0],p[1]};
          if (GridTest(exterior, pipoint)) {
            r_fp.add_point(p[0], p[1], 1, RasterTools::MAX);
          } else {
            r_fp.add_point(p[0], p[1], 0, RasterTools::MAX);
          }
          for (auto& hole : holes) {
          if (GridTest(hole, pipoint)) {
            r_fp.add_point(p[0], p[1], 1, RasterTools::MAX);
          } else {
            r_fp.add_point(p[0], p[1], 0, RasterTools::MAX);
          }
          }
          delete pipoint;
        }
      }
      delete exterior;
      for (auto& hole: holes) delete hole;
    }

    for(auto& p : points) {
      if (r_max.check_point(p[0], p[1])) {
        r_max.add_point(p[0], p[1], p[2], RasterTools::MAX);
        r_min.add_point(p[0], p[1], p[2], RasterTools::MIN);
        buckets[ r_max.getLinearCoord(p[0],p[1]) ].push_back(p[2]);
      }
    }
    
    PointCollection grid_points;
    vec1f values;
    double nodata = r_max.getNoDataVal();
    
    geoflow::Image I_max;
    I_max.dim_x = r_max.dimx_;
    I_max.dim_y = r_max.dimy_;
    I_max.min_x = r_max.minx_;
    I_max.min_y = r_max.miny_;
    I_max.cellsize = r_max.cellSize_;
    I_max.nodataval = r_max.noDataVal_;
    I_max.array = *r_max.vals_;
    geoflow::Image I_min(I_max);
    I_min.nodataval = r_min.noDataVal_;
    I_min.array = *r_min.vals_;
    geoflow::Image I_fp(I_max);
    I_fp.array = *r_fp.vals_;
    geoflow::Image I_cnt(I_max), I_med(I_max), I_avg(I_max), I_var(I_max);

    for(size_t i=0; i<r_max.dimy_ ; ++i) {
      for(size_t j=0; j<r_max.dimx_ ; ++j) {
        auto p = r_max.getPointFromRasterCoords(i,j);
        if (p[2]!=nodata) {
          grid_points.push_back(p);
          values.push_back(p[2]);
        }
        auto lc = r_max.getLinearCoord(i,j);
        auto& buck = buckets.at( lc );
        if (buck.size() == 0) {
          I_cnt.array[lc] = I_cnt.nodataval;
          I_med.array[lc] = I_med.nodataval;
          I_avg.array[lc] = I_avg.nodataval;
          I_var.array[lc] = I_var.nodataval;
        } else {
          std::sort(buck.begin(), buck.end());
          I_cnt.array[lc] = buck.size();
          I_med.array[lc] = buck[ buck.size()/2 ];
          I_avg.array[lc] = std::accumulate(buck.begin(), buck.end(), 0) / buck.size();
          int ssum = 0;
          for(auto& z : buck) {
            ssum += std::pow(z-I_avg.array[lc], 2);
          }
          I_var.array[lc] = ssum / buck.size();
        }
      }
    }
    
    output("image_fp").set(I_fp);
    output("image_max").set(I_max);
    output("image_min").set(I_min);
    output("image_cnt").set(I_cnt);
    output("image_med").set(I_med);
    output("image_avg").set(I_avg);
    output("image_var").set(I_var);
    output("heightfield").set(r_max);
    output("values").set(values);
    output("grid_points").set(grid_points);
  }


  void BuildingRasteriseNode::process() {
    auto& points = input("points").get<PointCollection&>();
    auto& ground_points = input("ground_points").get<PointCollection&>();
    auto h_ground = input("h_ground").get<float>();
    auto& footprint = input("footprint").get<LinearRing&>();

    auto box = footprint.box();
    auto boxmin = box.min();
    auto boxmax = box.max();

    RasterTools::Raster r(cellsize, boxmin[0]-cellsize, boxmax[0]+cellsize, boxmin[1]-cellsize, boxmax[1]+cellsize);
    r.prefill_arrays(RasterTools::MAX);

    // point in polygon to set ground plane outside footprint
    if (use_ground_pts) {
      auto exterior = build_grid(footprint);
      std::vector<pGridSet> holes;
      for (auto& hole : footprint.interior_rings()) {
        holes.push_back(build_grid(hole));
      }
      
      for (size_t col = 0; col < r.dimx_; ++col) {
        for (size_t row = 0; row < r.dimy_; ++row) {
          auto p = r.getPointFromRasterCoords(col, row);
          pPipoint pipoint = new Pipoint{p[0],p[1]};
          if (!GridTest(exterior, pipoint)) {
            if (r.add_point(p[0], p[1], h_ground, RasterTools::MAX)) {
              // ++data_pixel_cnt; //only count new cells (that were not written to before)
            }
          }
          for (auto& hole : holes) {
            if (GridTest(hole, pipoint)) {
              r.add_point(p[0], p[1], h_ground, RasterTools::MAX);
            }
          }
          delete pipoint;
        }
      }
      delete exterior;
      for (auto& hole: holes) delete hole;
    }

    if (use_tin) {
      PointCloud25DTriangulator pdt(
        do_area_thres,
        area_thres,
        do_angle_thres,
        angle_thres,
        do_normal_angle_thres,
        normal_angle_thres,
        do_len_thres,
        len_thres
      );
      for(auto& p : points) {
        pdt.insert_point(p);
      }
      if (use_ground_pts) {
        for(auto& p : ground_points) {
          pdt.insert_point(p);
        }
      }
      pdt.mark_big_triangles();
      pdt.remove_marked();
      pdt.mark_big_triangles();
      TriangleCollection triangles;
      pdt.get_triangles(triangles);
      
      float z_interpolate;
      for(size_t col=0; col<r.dimx_; ++col) {
        for(size_t row=0; row<r.dimy_; ++row) {
          // do linear TIN interpolation
          if (r.isNoData(col, row)) {
            auto p = r.getPointFromRasterCoords(col, row);
            if( pdt.get_z(p, z_interpolate) ) {
              r.set_val(col, row, z_interpolate);
            }
          }
        }
      }
    } else {
      for(auto& p : points) {
        if (r.check_point(p[0], p[1])) {
          r.add_point(p[0], p[1], p[2], RasterTools::MAX);
        }
      }
      if (use_ground_pts) {
        for(auto& p : ground_points) {
          if (r.check_point(p[0], p[1])) {
            if (r.isNoData(p[0], p[1])) {
              r.add_point(p[0], p[1], p[2], RasterTools::MAX);
            }
          }
        }
      }
    }

    PointCollection grid_points;
    grid_points.reserve(r.dimx_*r.dimy_);
    vec1f values;
    values.reserve(r.dimx_*r.dimy_);
    double nodata = r.getNoDataVal();
    for(size_t i=0; i<r.dimx_ ; ++i) {
      for(size_t j=0; j<r.dimy_ ; ++j) {
        auto p = r.getPointFromRasterCoords(i,j);
        if (p[2]!=nodata) {
          grid_points.emplace_back(p);
          values.emplace_back(p[2]);
        }
      }
    }
    geoflow::Image I;
    I.dim_x = r.dimx_;
    I.dim_y = r.dimy_;
    I.min_x = r.minx_;
    I.min_y = r.miny_;
    I.cellsize = r.cellSize_;
    I.nodataval = r.noDataVal_;
    I.array = *r.vals_;
    output("image").set(I);
    output("heightfield").set(r);
    output("values").set(values);
    output("grid_points").set(grid_points);
  }

  void RoofPartitionRasteriseNode::rasterise_arrangement(Arrangement_2& arr, RasterTools::Raster& r, size_t& data_pixel_cnt, bool use_planes) {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    for (auto face: arr.face_handles()) {
      if (face->data().in_footprint){ 
        
        auto& plane = face->data().plane;
        double z_interpolate = face->data().elevation_70p;
        LinearRing polygon;
        arrangementface_to_polygon(face, polygon);
        auto box = polygon.box();
        auto bb_min = box.min();
        auto bb_max = box.max();
        auto cr_min = r.getColRowCoord(bb_min[0], bb_min[1]);
        auto cr_max = r.getColRowCoord(bb_max[0], bb_max[1]);

        auto points_inside = r.rasterise_polygon(polygon, cr_min, cr_max);
        for (auto& p : points_inside) {
          if (use_planes)
            z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
          if (r.add_point(p[0], p[1], z_interpolate, RasterTools::MAX)) {
            ++data_pixel_cnt; //only count new cells (that were not written to before)
          }
        }
        // do plane projection
        // auto& plane = pts_per_roofplane[roofplane_ids[i]].first;
        // double z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
        // r.add_point(p[0], p[1], z_interpolate, RasterTools::MAX);
      }
    }
  }
  void RoofPartitionRasteriseNode::process(){
    auto arr = input("arrangement").get<Arrangement_2>();
    auto h_ground = input("h_ground").get<float>();

    auto unbounded_face = arr.unbounded_face();
    auto& footprint = input("footprint").get<LinearRing&>();

    auto box = footprint.box();
    auto boxmin = box.min();
    auto boxmax = box.max();
    RasterTools::Raster r = RasterTools::Raster(cellsize, boxmin[0]-0.5, boxmax[0]+0.5, boxmin[1]-0.5, boxmax[1]+0.5);
    r.prefill_arrays(RasterTools::MAX);
    size_t roofdata_area_cnt = 0;
    rasterise_arrangement(arr, r, roofdata_area_cnt, use_planes_);

    float volume = 0;
    float cellarea = cellsize*cellsize;
    for (size_t col = 0; col < r.dimx_; ++col) {
      for (size_t row = 0; row < r.dimy_; ++row) {
        if (!r.isNoData(col, row)) {
          volume += cellarea * std::fabs(r.get_val(col, row) - h_ground);
        }
      }
    }
    output("volume").set(volume);

    geoflow::Image I;
    I.dim_x = r.dimx_;
    I.dim_y = r.dimy_;
    I.min_x = r.minx_;
    I.min_y = r.miny_;
    I.cellsize = r.cellSize_;
    I.nodataval = r.noDataVal_;
    I.array = *r.vals_;
    output("image").set(I);
  }

  void rasterise_ring(LinearRing& polygon, RasterTools::Raster& r) {
    typedef double                      FT;
    typedef CGAL::Simple_cartesian<FT>  K;

    K::Plane_3 plane;
    std::vector<K::Point_3> pts;
    for (auto& p : polygon) {
      pts.push_back(K::Point_3(p[0], p[1], p[2]));
    }
    // fit plane to polygon pts
    linear_least_squares_fitting_3(pts.begin(),pts.end(),plane,CGAL::Dimension_tag<0>());

    auto box = polygon.box();
    auto bb_min = box.min();
    auto bb_max = box.max();
    auto cr_min = r.getColRowCoord(bb_min[0], bb_min[1]);
    auto cr_max = r.getColRowCoord(bb_max[0], bb_max[1]);

    auto points_inside = r.rasterise_polygon(polygon, cr_min, cr_max);
    for (auto& p : points_inside) {
      float z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
      r.add_point(p[0], p[1], z_interpolate, RasterTools::MAX);
    }

  }
  
  void calculate_h_attr(gfSingleFeatureInputTerminal& roofparts, gfMultiFeatureOutputTerminal& h_attr, RasterTools::Raster& r_lod22, float z_offset) {
    h_attr.add_vector("h_min", typeid(float));
    h_attr.add_vector("h_max", typeid(float));
    h_attr.add_vector("h_50p", typeid(float));
    h_attr.add_vector("h_70p", typeid(float));
    for (size_t i=0; i< roofparts.size(); ++i) {
      auto polygon = roofparts.get<LinearRing>(i);
      auto box = polygon.box();
      auto bb_min = box.min();
      auto bb_max = box.max();
      auto cr_min = r_lod22.getColRowCoord(bb_min[0], bb_min[1]);
      auto cr_max = r_lod22.getColRowCoord(bb_max[0], bb_max[1]);
      auto part_points = r_lod22.rasterise_polygon(polygon, cr_min, cr_max, false);

      if(part_points.size()==0) {
        part_points.insert(part_points.begin(), polygon.begin(), polygon.end());
      }
      std::sort(part_points.begin(), part_points.end(), [](auto& p1, auto& p2) {
        return p1[2] < p2[2];
      });

      size_t datasize = part_points.size();
      int elevation_id = std::floor(0.5*float(datasize-1));
      h_attr.sub_terminal("h_50p").push_back(part_points[elevation_id][2] + z_offset);
      elevation_id = std::floor(0.7*float(datasize-1));
      h_attr.sub_terminal("h_70p").push_back(part_points[elevation_id][2] + z_offset);
      h_attr.sub_terminal("h_min").push_back(part_points[0][2] + z_offset);
      h_attr.sub_terminal("h_max").push_back(part_points[datasize-1][2] + z_offset);
    }
  }

  typedef CGAL::Simple_cartesian<float> CF;
  CF::Vector_3 calculate_normal_cf(const LinearRing& ring)
  {
    float x=0, y=0, z=0;
    for (size_t i = 0; i < ring.size(); ++i) {
      const auto &curr = ring[i];
      const auto &next = ring[(i + 1) % ring.size()];
      x += (curr[1] - next[1]) * (curr[2] + next[2]);
      y += (curr[2] - next[2]) * (curr[0] + next[0]);
      z += (curr[0] - next[0]) * (curr[1] + next[1]);
    }
    CF::Vector_3 n(x, y, z);
    return n / CGAL::approximate_sqrt(n.squared_length());
  }
  
  static constexpr double pi = 3.14159265358979323846;
  static const CF::Vector_3 up = CF::Vector_3(0,0,1);

  void RoofPartition3DBAGRasteriseNode::process(){
    auto& lod12_roofparts = input("lod12_roofparts");//.get<LinearRing>();
    auto& lod13_roofparts = input("lod13_roofparts");//.get<LinearRing>();
    auto& lod22_roofparts = input("lod22_roofparts");//.get<LinearRing>();

    auto& lod12_hattr = poly_output("lod12_hattr");
    auto& lod13_hattr = poly_output("lod13_hattr");
    auto& lod22_hattr = poly_output("lod22_hattr");

    Box box;
    for (size_t i=0; i< lod22_roofparts.size(); ++i) {
      auto rpart = lod22_roofparts.get<LinearRing>(i);
      box.add(rpart.box());
    }
    auto boxmin = box.min();
    auto boxmax = box.max();

    RasterTools::Raster r_lod22 = RasterTools::Raster(cellsize, boxmin[0]-0.5, boxmax[0]+0.5, boxmin[1]-0.5, boxmax[1]+0.5);
    r_lod22.prefill_arrays(RasterTools::MAX);
    for (size_t i=0; i< lod22_roofparts.size(); ++i) {
      auto ring = lod22_roofparts.get<LinearRing>(i);
      rasterise_ring(ring, r_lod22);
    }
    auto z_offset = (*manager.data_offset())[2];
    calculate_h_attr(lod12_roofparts, lod12_hattr, r_lod22, z_offset);
    calculate_h_attr(lod13_roofparts, lod13_hattr, r_lod22, z_offset);
    calculate_h_attr(lod22_roofparts, lod22_hattr, r_lod22, z_offset);

    // compute inclination and azimuth. Both in degrees.
    lod22_hattr.add_vector("b3_hellingshoek", typeid(float));
    lod22_hattr.add_vector("b3_azimut", typeid(float));
    for (size_t i=0; i<lod22_roofparts.size(); ++i) {
      auto ring = lod22_roofparts.get<LinearRing>(i);
      auto n = calculate_normal_cf(ring);
      float azimuth, slope;
      if (std::isnan(n.x()) || std::isnan(n.y()) || std::isnan(n.z())){
        azimuth = std::numeric_limits<float>::quiet_NaN();
        slope = std::numeric_limits<float>::quiet_NaN();
      } else {
        slope = CGAL::to_double(CGAL::approximate_angle(n, up));
        
        // calculate azimuth from arctan2 (https://en.cppreference.com/w/cpp/numeric/math/atan2)
        // ie. subtract pi/2, multiply by -1 and then add 2 pi if result is negative (4th quadrant)
        azimuth = -1 * ( std::atan2(n.y(), n.x()) - pi/2 );
        if (azimuth<0) {
          azimuth = 2*pi + azimuth;
        }

        // convert to degrees
        azimuth = azimuth * (180/pi);
      }
      // push attributes
      lod22_hattr.sub_terminal("b3_hellingshoek").push_back(slope);
      lod22_hattr.sub_terminal("b3_azimut").push_back(azimuth);
    }
  }

  void SegmentRasteriseNode::rasterise_input(gfSingleFeatureInputTerminal& input_triangles, RasterTools::Raster& r, size_t& data_pixel_cnt) {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    for(size_t i=0; i< input_triangles.size(); ++i) {
      // auto polygon = alpha_rings.get<LinearRing>(i);
      auto& triangle_collection = input_triangles.get<TriangleCollection>(i);
      // auto dt = alpha_dts[i];
      // as::Triangulation_2 T;
      // auto points = pts_per_roofplane[roofplane_ids[i]].second;
      // auto plane = pts_per_roofplane[roofplane_ids[i]].first;
      // T.insert(points.begin(), points.end());
      // auto as = as::Alpha_shape_2(T, as::FT(thres_alpha), as::Alpha_shape_2::GENERAL);
      for (auto& triangle : triangle_collection) {
        CGAL::Plane_3<K> plane(
          K::Point_3(triangle[0][0], triangle[0][1], triangle[0][2]),
          K::Point_3(triangle[1][0], triangle[1][1], triangle[1][2]),
          K::Point_3(triangle[2][0], triangle[2][1], triangle[2][2])
        );
        geoflow::Box box;
        for (auto& p : triangle) {
          box.add(p);
        }
        auto bb_min = box.min();
        auto bb_max = box.max();
        auto cr_min = r.getColRowCoord(bb_min[0], bb_min[1]);
        auto cr_max = r.getColRowCoord(bb_max[0], bb_max[1]);

        auto points_inside = r.rasterise_polygon(triangle, cr_min, cr_max);
        for (auto& p : points_inside) {
          double z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
          if (r.add_point(p[0], p[1], z_interpolate, RasterTools::MAX)) {
            ++data_pixel_cnt; //only count new cells (that were not written to before)
          }
        }
        // do plane projection
        // auto& plane = pts_per_roofplane[roofplane_ids[i]].first;
        // double z_interpolate = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
        // r.add_point(p[0], p[1], z_interpolate, RasterTools::MAX);
      }
    }
  }
  void SegmentRasteriseNode::process() {
    // auto& alpha_rings = vector_input("alpha_rings");
    auto& ground_triangles = vector_input("ground_triangles");
    auto& triangles = vector_input("triangles");
    // auto& alpha_dts = input("alpha_dts").get<std::vector<as::Triangulation_2>&>();
    // auto& roofplane_ids = input("roofplane_ids").get<vec1i&>();
    // auto& pts_per_roofplane = input("pts_per_roofplane").get<IndexedPlanesWithPoints&>();
    RasterTools::Raster r;

    Box box;
    for(size_t i=0; i< triangles.size(); ++i){
      auto triangle_collection = triangles.get<TriangleCollection>(i);
      box.add(triangle_collection.box());
    }
    if(use_ground) {
      for(size_t i=0; i< ground_triangles.size(); ++i){
        auto triangle_collection = ground_triangles.get<TriangleCollection>(i);
        box.add(triangle_collection.box());
      }
    }
    auto boxmin = box.min();
    auto boxmax = box.max();
    auto pixel_limit = megapixel_limit * 1E6;
    while(true) {
      auto dimx = (boxmax[0]-boxmin[0])/cellsize + 1;
      auto dimy = (boxmax[1]-boxmin[1])/cellsize + 1;
      if(dimx*dimy > pixel_limit) {
        cellsize*=2;
      } else break;
    }
    r = RasterTools::Raster(cellsize, boxmin[0]-0.5, boxmax[0]+0.5, boxmin[1]-0.5, boxmax[1]+0.5);
    r.prefill_arrays(RasterTools::MAX);

    size_t roofdata_area_cnt = 0, grounddata_area_cnt = 0;
    rasterise_input(triangles, r, roofdata_area_cnt);
    if(use_ground)
      rasterise_input(ground_triangles, r, grounddata_area_cnt);

    if (fill_nodata_) r.fill_nn(fill_nodata_window_size_);

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
    output("data_area").set(float(roofdata_area_cnt)*cellsize*cellsize);
    output("heightfield").set(r);
    output("values").set(values);
    output("grid_points").set(grid_points);
  }

  void RasterMergerNode::process() {
    auto rasterbase = input("rasterbase").get<RasterTools::Raster>();
    auto& rasteradd = input("rasteradd").get<RasterTools::Raster>();
    
    for(size_t i=0; i<rasteradd.dimx_ ; ++i) {
      for(size_t j=0; j<rasteradd.dimy_ ; ++j) {
        auto p = rasteradd.getPointFromRasterCoords(i,j);
        rasterbase.add_point(p[0], p[1], p[2], RasterTools::MAX);
      }
    }

    output("raster").set(rasterbase);
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