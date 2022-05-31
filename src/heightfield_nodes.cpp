// triangulation
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Projection_traits_xy_3.h>

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
          box.add({float(p.x()), float(p.y()), float(p.z())});
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

    RasterTools::Raster r(cellsize, boxmin[0], boxmax[0], boxmin[1], boxmax[1]);
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
        r.add_point(p[0], p[1], p[2], RasterTools::MAX);
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

    auto unbounded_face = arr.unbounded_face();
    auto& footprint = input("footprint").get<LinearRing&>();

    auto box = footprint.box();
    auto boxmin = box.min();
    auto boxmax = box.max();
    RasterTools::Raster r = RasterTools::Raster(cellsize, boxmin[0]-0.5, boxmax[0]+0.5, boxmin[1]-0.5, boxmax[1]+0.5);
    r.prefill_arrays(RasterTools::MAX);
    size_t roofdata_area_cnt = 0;
    rasterise_arrangement(arr, r, roofdata_area_cnt, use_planes_);

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