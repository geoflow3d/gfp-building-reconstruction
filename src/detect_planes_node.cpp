#include "stepedge_nodes.hpp"
// #include <CGAL/Regularization/regularize_planes.h>

#include "plane_detect.hpp"

namespace geoflow::nodes::stepedge {

  void DetectPlanesNode::process() {
    auto points = input("points").get<PointCollection>();

    // convert to cgal points with attributes
    PNL_vector pnl_points;
    for (auto& p : points) {
      PNL pv;
      boost::get<0>(pv) = Point(p[0], p[1], p[2]);
      boost::get<2>(pv) = 0;
      boost::get<3>(pv) = 0;
      boost::get<9>(pv) = 0;
      pnl_points.push_back(pv);
    }
    // estimate normals
    if (points.size()) {
      CGAL::pca_estimate_normals<Concurrency_tag>(
        pnl_points, metrics_normal_k,
        CGAL::parameters::point_map(Point_map()).
        normal_map(Normal_map())
      );
    }
    // orient normals upwards
    auto up = Vector(0,0,1);
    for ( auto& pv : pnl_points) {
      auto &n = boost::get<1>(pv);
      if (n*up<0) 
        boost::get<1>(pv) = -n;
    }

    // convert to lists required by the planedetector class
    // size_t i=0;
    PointCollection points_vec;
    vec3f normals_vec;
    points_vec.reserve(points.size());
    for (auto &pt : pnl_points) {
      auto& p = boost::get<0>(pt);
      auto& n = boost::get<1>(pt);
      points_vec.push_back(
        {float(CGAL::to_double(p.x())), float(CGAL::to_double(p.y())), float(CGAL::to_double(p.z()))}
      );
      normals_vec.push_back(
        {float(CGAL::to_double(n.x())), float(CGAL::to_double(n.y())), float(CGAL::to_double(n.z()))}
      );
    }
    // perform plane detection
    planedect::PlaneDS PDS(points_vec, normals_vec, metrics_plane_k);
    planedect::DistAndNormalTester DNTester(
      metrics_plane_epsilon * metrics_plane_epsilon,
      metrics_plane_normal_threshold,
      n_refit
    );
    regiongrower::RegionGrower<planedect::PlaneDS, planedect::PlaneRegion> R;
    R.min_segment_count = metrics_plane_min_points;
    if(points.size()>metrics_plane_min_points)
      R.grow_regions(PDS, DNTester);

    
    // if (regularise_planes) {
    //   // Regularize detected planes.
    //   CGAL::regularize_planes(points,
    //                           Point_map(),
    //                           planes,
    //                           CGAL::Shape_detection::Plane_map<Traits>(),
    //                           CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
    //                           true,  // regularize parallelism
    //                           true,  // regularize orthogonality
    //                           false, // do not regularize coplanarity
    //                           true,  // regularize Z-symmetry (default)
    // }


    // classify horizontal/vertical planes using plane normals
    IndexedPlanesWithPoints pts_per_roofplane;
    size_t horiz_roofplane_cnt=0;
    size_t slant_roofplane_cnt=0;
    if (only_horizontal) pts_per_roofplane[-1].second = std::vector<Point>();
    size_t horiz_pt_cnt=0, total_pt_cnt=0, wall_pt_cnt=0, unsegmented_pt_cnt=0;
    vec1f roof_elevations;
    
    for(auto region: R.regions){
      auto& plane = region.plane;
      output("planes").push_back(plane);
      Vector n = plane.orthogonal_vector();
      // this dot product is close to 0 for vertical planes
      auto horizontality = CGAL::abs(n*Vector(0,0,1));
      bool is_wall = horizontality < metrics_is_wall_threshold;
      bool is_horizontal = horizontality > metrics_is_horizontal_threshold;
      // put slanted surface points at index -1 if we care only about horzontal surfaces
      if (!is_wall) {
        std::vector<Point> segpts;
        for (auto& i : region.inliers) {
          segpts.push_back(boost::get<0>(pnl_points[i]));
          roof_elevations.push_back(float(boost::get<0>(pnl_points[i]).z()));
        }
        total_pt_cnt += segpts.size();
        if (!only_horizontal ||
            (only_horizontal && is_horizontal)) {
          pts_per_roofplane[region.get_region_id()].second = segpts;
          pts_per_roofplane[region.get_region_id()].first = plane;
        } else if (!is_horizontal) {
          pts_per_roofplane[-1].second.insert(
            pts_per_roofplane[-1].second.end(),
            segpts.begin(),
            segpts.end()
          );
        } 
        if (is_horizontal) {
          horiz_pt_cnt += segpts.size();
        }
      } else { // is_wall
        wall_pt_cnt = region.inliers.size();
      }
      if (is_horizontal)
        ++horiz_roofplane_cnt;
      else if (!is_wall && !is_horizontal)
        ++slant_roofplane_cnt;

      for (size_t& i : region.inliers) {
        boost::get<2>(pnl_points[i]) = region.get_region_id();
        boost::get<3>(pnl_points[i]) = is_wall;
        boost::get<9>(pnl_points[i]) = is_horizontal;
      }
    }

    bool b_is_horizontal = float(horiz_pt_cnt)/float(total_pt_cnt) > horiz_min_count;
    // int roof_type=-2; // as built: -2=undefined; -1=no pts; 0=LOD1, 1=LOD1.3, 2=LOD2
    std::string roof_type = "no planes";
    if (R.regions.size()==0) {
      // roof_type=-1;
      roof_type = "no points";
    } else if (horiz_roofplane_cnt==1 && slant_roofplane_cnt==0){
      // roof_type=0;
      roof_type = "horizontal";
    } else if (b_is_horizontal){
      // roof_type=1;
      roof_type = "multiple horizontal";
    } else if (slant_roofplane_cnt > 0) {
      // roof_type=2;
      roof_type = "slanted";
    }

    if (roof_elevations.size()) {
      output("roof_elevation_70p").set(compute_percentile(roof_elevations, 0.7));
      output("roof_elevation_50p").set(compute_percentile(roof_elevations, 0.5));
      output("roof_elevation_min").set(compute_percentile(roof_elevations, 0.0));
      output("roof_elevation_max").set(compute_percentile(roof_elevations, 1.0));
    } else {
      output("roof_elevation_70p").set_from_any(std::any());
      output("roof_elevation_50p").set_from_any(std::any());
      output("roof_elevation_min").set_from_any(std::any());
      output("roof_elevation_max").set_from_any(std::any());
    }

    output("roof_type").set(roof_type);
    output("horiz_roofplane_cnt").set(float(horiz_roofplane_cnt));
    output("slant_roofplane_cnt").set(float(slant_roofplane_cnt));
    output("total_roofplane_cnt").set(int(horiz_roofplane_cnt+slant_roofplane_cnt));
    output("roof_pt_cnt").set((int)total_pt_cnt);
    output("wall_pt_cnt").set((int)wall_pt_cnt);

    vec1i plane_id, is_wall, is_horizontal;
    for(auto& p : pnl_points) {
      auto pid = boost::get<2>(p);
      if (pid==0) ++unsegmented_pt_cnt;
      plane_id.push_back(pid);
      is_wall.push_back(boost::get<3>(p));
      is_horizontal.push_back(boost::get<9>(p));
    }
    output("unsegmented_pt_cnt").set((int)unsegmented_pt_cnt);
    output("pts_per_roofplane").set(pts_per_roofplane);
    output("plane_id").set(plane_id);
    output("is_wall").set(is_wall);
    output("is_horizontal").set(is_horizontal);

    output("plane_adj").set(R.adjacencies);
  }

}