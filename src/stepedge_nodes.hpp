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
#pragma once

#include <geoflow/common.hpp>
#include <geoflow/geoflow.hpp>

#include "arrangement.hpp"
#include "point_edge.h"
#include "line_regulariser.hpp"
#include "Raster.h"

namespace geoflow::nodes::stepedge {

  float compute_percentile(std::vector<float>& z_vec, float percentile);

  typedef std::unordered_map<int, std::pair<Plane, std::vector<Point>>> IndexedPlanesWithPoints;

  typedef std::vector<Mesh> MultiMesh;

  class AlphaShapeNode:public Node {
    float thres_alpha = 0.25;
    bool extract_polygons = false;
    bool optimal_alpha = true;
    bool optimal_only_if_needed = true;
    public:
    using Node::Node;
    void init() override {
      // add_input("points", TT_any);
      add_input("pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
      add_input("skip", typeid(bool));
      add_vector_output("alpha_rings", typeid(LinearRing));
      add_output("edge_points", typeid(PointCollection));
      add_output("alpha_edges", typeid(LineStringCollection));
      add_vector_output("alpha_triangles", typeid(TriangleCollection));
      // add_output("alpha_dts", typeid(std::vector<as::Triangulation_2>));
      add_output("segment_ids", typeid(vec1i));
      add_output("boundary_points", typeid(PointCollection));
      add_output("roofplane_ids", typeid(vec1i));

      add_param(ParamFloat(thres_alpha, "thres_alpha", "thres_alpha"));
      add_param(ParamBool(optimal_alpha, "optimal_alpha", "optimal_alpha"));
      add_param(ParamBool(optimal_only_if_needed, "optimal_only_if_needed", "optimal_only_if_needed"));
      add_param(ParamBool(extract_polygons, "extract_polygons", "extract_polygons"));
    }
    
    bool inputs_valid() override {
      auto& skipInput = input("skip");
      if (skipInput.has_data() && input("pts_per_roofplane").has_data()) {
        return !skipInput.get<bool>();
      }
      return false;
    }

    // void before_gui(){
    //   auto param = std::get<ParamBool>(parameters.at("optimal_only_if_needed"));
    //   param.set_visible(optimal_alpha);
    // }
    void process() override;
  };

  class Ring2SegmentsNode:public Node {
    public:
    using Node::Node;
    void init() override {
      add_input("rings", typeid(LinearRingCollection));
      add_output("edge_segments", typeid(SegmentCollection));
      add_output("ring_idx", typeid(std::vector<std::vector<size_t>>));
    }
    void process() override;
  };

  class PointCloudMeanZNode:public Node {

    public:
    using Node::Node;
    void init() override {
      add_input("point_clouds", typeid(std::vector<PointCollection>));
      add_output("height", typeid(vec1f));
    }

    void process() override;
  };

  class PolygonExtruderNode:public Node {

    public:
    using Node::Node;
    void init() override {
      add_input("polygons", typeid(LinearRingCollection));
      add_input("heights", typeid(vec1f));
      add_output("rings_3d", typeid(LinearRingCollection));
      add_output("ring_types", typeid(vec1i));
    }

    void process() override;
  };

  class LOD1ExtruderNode:public Node {

    public:
    using Node::Node;
    void init() override {
      add_input("polygon", typeid(LinearRing));
      add_input("floor_elevation", typeid(float));
      add_input("roof_elevation", typeid(float));
      add_vector_output("3d_polygons", typeid(LinearRing));
      add_output("surface_types", typeid(vec1i));
      add_output("mesh", typeid(std::unordered_map<int, Mesh>));
    }

    // bool inputs_valid() override {
    //   return input("polygon").has_data() && input("floor_elevation").has_data();
    // }

    void process() override;
  };

  class Arr2LinearRingsNode:public Node {
    // bool invalid_rooftype = false;

    bool only_in_footprint = true;
    bool output_groundparts = false;
    int plane_id = 0;
    public:
    using Node::Node;
    void init() override {
      add_vector_input("arrangement", typeid(Arrangement_2));
      add_vector_input("groundparts", typeid(LinearRing));
      // add_vector_input("floor_elevation", typeid(float));
      // add_vector_input("mesh_error", typeid(float));
      // add_vector_input("roof_type", typeid(int));
      // add_vector_input("arr_complexity", typeid(int));
      add_poly_input("attributes", {typeid(bool), typeid(int), typeid(float), typeid(std::string), typeid(std::string), typeid(Date), typeid(Time), typeid(DateTime)});
      add_poly_output("attributes", {typeid(bool), typeid(int), typeid(float), typeid(std::string), typeid(std::string), typeid(Date), typeid(Time), typeid(DateTime)});
      add_vector_output("linear_rings", typeid(LinearRing));
      // add_vector_output("plane_a", typeid(float));
      // add_vector_output("plane_b", typeid(float));
      // add_vector_output("plane_c", typeid(float));
      // add_vector_output("plane_d", typeid(float));

      add_param(ParamBool(output_groundparts, "output_groundparts", "Also output the ground parts"));
    }
    bool inputs_valid() override {
      return vector_input("arrangement").has_data() && poly_input("attributes").has_data() && vector_input("arrangement").is_touched();
    }
    void process() override;
  };

  class Arr2LinearRingsDebugNode:public Node {
    public:
    using Node::Node;
    void init() override {
      add_vector_input("arrangement", typeid(Arrangement_2));
      add_poly_output("attributes", {typeid(bool), typeid(int), typeid(float), typeid(std::string)});
      add_vector_output("linear_rings", typeid(LinearRing));
      add_output("isolated_vertices", typeid(PointCollection));
    }
    void process() override;
  };

  class ArrExtruderNode:public Node {
    bool do_walls=true, do_roofs=true, do_floor=true;
    bool LoD2 = true;
    bool lod1_extrude_to_max_ = false;
    // float base_elevation = 0;
    float nodata_elevation = 3;
    int snap_tolerance_exp = 4;

    public:
    using Node::Node;
    void init() override {
      add_input("arrangement", typeid(Arrangement_2));
      add_input("floor_elevation", typeid(float));

      // add_output("normals_vec3f", typeid(vec3f), true);
      add_vector_output("labels", typeid(int)); // 0==ground, 1==roof, 2==outerwall, 3==innerwall
      add_vector_output("faces", typeid(LinearRing));
      add_vector_output("mesh", typeid(Mesh));
      add_output("multisolid", typeid(std::unordered_map<int, Mesh>));

      add_param(ParamBool(do_walls, "do_walls", "Do walls"));
      add_param(ParamBool(do_roofs, "do_roofs", "Do roofs"));
      add_param(ParamBool(do_floor, "do_floor", "Do floor"));
      add_param(ParamBool(LoD2, "LoD2", "LoD2 (uncheck for LoD1.3)"));
      add_param(ParamBool(lod1_extrude_to_max_, "lod1_extrude_to_max", "Extrude LoD1.x roofparts to 97th percentile height. If false the 70th percentile is used."));
      // ", ParamFloat( add_param("b,base_elevation, "Base elevation"));
      // add_param(ParamFloat(nodata_elevation, "nodata_elevation", "Nodata elevation"));
      add_param(ParamInt(snap_tolerance_exp, "snap_tolerance_exp", "Snap tolerance"));
    }
    void process() override;
  };

  // class ExtruderNode:public Node {
  //   bool do_walls=true, do_roofs=true;
  //   bool in_footprint = false;
  //   bool LoD2 = false;
  //   float base_elevation = 0;
  //   float nodata_elevation = 3;

  //   public:
  //   using Node::Node;
  //   void init() override {
  //     add_input("arrangement", typeid(Arrangement_2));
  //     add_output("cell_id_vec1i", typeid(vec1i));
  //     add_output("plane_id", typeid(vec1i));
  //     add_output("rms_errors", typeid(vec1f));
  //     add_output("max_errors", typeid(vec1f));
  //     add_output("elevations", typeid(vec1f));
  //     add_output("segment_coverages", typeid(vec1f));
  //     add_output("triangles", typeid(TriangleCollection));
  //     add_output("normals_vec3f", typeid(vec3f));
  //     add_output("labels_vec1i", typeid(vec1i)); // 0==ground, 1==roof, 2==outerwall, 3==innerwall
  //     add_output("face_ids", typeid(vec1i)); // 0==ground, 1==roof, 2==outerwall, 3==innerwall

  //     add_param(ParamBool(do_walls, "do_walls", "Do walls"));
  //     add_param(ParamBool(do_roofs, "do_roofs", "Do roofs"));
  //     add_param(ParamBool(in_footprint, "in_footprint", "in_footprint"));
  //     add_param(ParamBool(LoD2, "LoD2", "LoD2"));
  //     add_param(ParamFloat(base_elevation, "base_elevation", "Base elevation"));
  //     add_param(ParamFloat(nodata_elevation, "nodata_elevation", "Nodata elevation"));
  //   }
  //   void process() override;
  // };

  class PolygonGrowerNode:public Node {
    float extension = 0.1;
    public:
    using Node::Node;
    void init() override {
      add_input("rings", typeid(LinearRingCollection));
      add_output("rings", typeid(LinearRingCollection));

      add_param(ParamFloat(extension, "extension", "Extension length"));
    }
    void process() override;
  };

  // class BuildArrFromRingsExactNode:public Node {
  //   bool extrude_unsegmented = true;
  //   float extrude_mindensity = 5;
  //   float z_percentile = 0.;
  //   float rel_area_thres = 0.1;
  //   bool flood_to_unsegmented = true;
  //   bool dissolve_edges = true;
  //   bool dissolve_stepedges = true;
  //   float step_height_threshold = 1.0;
  //   bool snap_clean = true;
  //   bool snap_clean_fp = false;
  //   bool snap_detect_only = false;
  //   float snap_dist = 1.0;

  //   public:
  //   bool arr_is_valid=false;
  //   int vcount, ecount;

  //   using Node::Node;
  //   void init() override {
  //     add_input("rings", typeid(std::unordered_map<size_t, linereg::Polygon_2>));
  //     add_input("pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
  //     add_input("footprint", {typeid(linereg::Polygon_2), typeid(LinearRing)});
  //     add_output("noseg_area_a", typeid(float));
  //     add_output("noseg_area_r", typeid(float));
  //     add_output("arrangement", typeid(Arrangement_2));
  //     add_output("arr_segments", typeid(LineStringCollection));
  //     add_output("snap_to_e", typeid(SegmentCollection));
  //     add_output("snap_to_v", typeid(PointCollection));
  //     add_output("snap_v", typeid(PointCollection));
  //     add_output("snap_fp_to_e", typeid(SegmentCollection));
  //     add_output("snap_fp_to_v", typeid(PointCollection));
  //     add_output("snap_fp_v", typeid(PointCollection));

  //     add_param(ParamBool(extrude_unsegmented, "extrude_unsegmented", "extrude_unsegmented"));
  //     add_param(ParamBoundedFloat(extrude_mindensity, 1, 20, "extrude_mindensity",  "Extrude min density"));
  //     add_param(ParamBoundedFloat(z_percentile, 0, 1, "z_percentile",  "Elevation percentile"));
  //     add_param(ParamBoundedFloat(rel_area_thres, 0.01, 1,  "rel_area_thres", "Preserve split ring area"));
  //     add_param(ParamBool(flood_to_unsegmented, "flood_to_unsegmented", "flood_to_unsegmented"));
  //     add_param(ParamBool(dissolve_edges, "dissolve_edges", "dissolve_edges"));
  //     add_param(ParamBool(dissolve_stepedges, "dissolve_stepedges", "dissolve_stepedges"));
  //     add_param(ParamBoundedFloat(step_height_threshold, 0, 10, "step_height_threshold",  "step_height_threshold"));
  //     add_param(ParamBool(snap_clean, "snap_clean", "Snap"));
  //     add_param(ParamBool(snap_clean_fp, "snap_clean_fp", "Snap fp"));
  //     add_param(ParamBool(snap_detect_only, "snap_detect_only", "snap_detect_only"));
  //     add_param(ParamBoundedFloat(snap_dist, 0.01, 5, "snap_dist", "Snap distance"));
  //   }
  //   // std::string info() {
  //   //   ImGui::Text("Arrangement valid? %s", arr_is_valid? "yes" : "no");
  //   //   ImGui::Text("vcount: %d, ecount: %d", vcount, ecount);
  //   // }
  //   // void arr_snapclean(Arrangement_2& arr);
  //   void arr_snapclean_from_fp(Arrangement_2& arr);
  //   void arr_process(Arrangement_2& arr);
  //   void arr_assign_pts_to_unsegmented(Arrangement_2& arr, std::vector<Point>& points);
  //   void process() override;
  // };

  class BuildArrFromLinesNode:public Node {
    // float rel_area_thres = 0.1;
    int max_arr_complexity = 400;
    int dist_threshold_exp = 4;
    float fp_extension = 0.0;
    bool insert_with_snap = false;
    bool insert_lines = true;
    // int angle_threshold_exp = 5;
    // bool snap_clean = true;
    // bool snap_detect_only = false;
    // float snap_dist = 1.0;
    public:

    using Node::Node;
    void init() override {
      add_vector_input("lines", {typeid(Segment), typeid(linereg::Segment_2)});
      add_input("footprint", {typeid(linereg::Polygon_with_holes_2), typeid(LinearRing)});
      add_output("arrangement", typeid(Arrangement_2));
      add_output("arr_complexity", typeid(int));

      // add_param(ParamBoundedFloat(rel_area_thres, 0.01, 1,  "rel_area_thres", "Preserve split ring area"));
      add_param(ParamBoundedFloat(fp_extension, 0.0, 0.01,  "fp_extension", "extend each footprint segment on both sides with this distance. NB can interfere with in_footprint detection!"));
      add_param(ParamInt(max_arr_complexity, "max_arr_complexity", "Maximum nr of lines"));
      add_param(ParamBoundedInt(dist_threshold_exp, 0, 15, "dist_threshold_exp", "10-base exponent to set distance threshold. Eg a value of 2 yields a value of 10^(-2) = 0.01"));
      // add_param(ParamBoundedInt(angle_threshold_exp, 0, 15, "angle_threshold_exp", "10-base exponent to set angle threshold in degrees. Eg a value of 2 yields a value of 10^(-2) = 0.01"));

      add_param(ParamBool(insert_with_snap, "insert_with_snap", "Use custom insert function that aims to prevent duplicate points etc"));
      add_param(ParamBool(insert_lines, "insert_lines", "Insert lines from the lines input terminal"));
      // add_param(ParamBool(snap_clean, "snap_clean", "Snap"));
      // add_param(ParamBool(snap_detect_only, "snap_detect_only", "snap_detect_only"));
      // add_param(ParamBoundedFloat(snap_dist, 0.01, 5, "snap_dist", "Snap distance"));
    }
    bool inputs_valid() override {
      for (auto& [name,iT] : input_terminals) {
        if (name == "lines" && !insert_lines)
          continue;
        else if (!iT->is_touched())
          return false;
      }
        
      return true;
    }
    void process() override;
  };

  class OptimiseArrangmentGridNode:public Node {
    float data_multiplier = 50.0;
    float smoothness_multiplier = 1.0;
    bool preset_labels = false;
    bool do_normalise = true;
    int n_iterations = 3;
    int graph_cut_impl = 2;
    bool use_ground = true;
    bool label_ground_outside_fp = false;

    float z_percentile = 0.9;

    public:
    using Node::Node;
    void init() override {
      add_input("arrangement", typeid(Arrangement_2));
      add_input("heightfield", typeid(RasterTools::Raster));
      add_input("pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
      add_input("ground_pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
      add_output("arrangement", typeid(Arrangement_2));
      add_vector_output("groundparts", typeid(LinearRing));

      add_param(ParamBoundedInt(graph_cut_impl, 0, 2, "graph_cut_impl", "Graph cut implementation"));
      add_param(ParamBoundedInt(n_iterations, 0, 100, "n_iterations", "Number of iterations"));
      add_param(ParamBoundedFloat(data_multiplier, 0.001, 100,  "data_multiplier", "Multiplier on data term"));
      add_param(ParamBoundedFloat(smoothness_multiplier, 0.001, 100,  "smoothness_multiplier", "Multiplier on smoothness term"));
      add_param(ParamBool(preset_labels, "preset_labels", "Preset face labels"));
      add_param(ParamBool(do_normalise, "do_normalise", "Normalise weights"));
      add_param(ParamBoundedFloat(z_percentile, 0, 1, "z_percentile",  "Elevation percentile"));
      add_param(ParamBool(use_ground, "use_ground", "Use ground planes in optimisation"));
      add_param(ParamBool(label_ground_outside_fp, "label_ground_outside_fp", "Label faces that get assigned a ground plane as begin outside the footprint"));
    }
    bool inputs_valid() override {
      for (auto& [name,iT] : input_terminals) {
        if (name == "ground_pts_per_roofplane" && !use_ground) {
          continue;
        } else if (!iT->is_touched()) {
          return false;
        }
      }
      return true;
    }
    void process() override;
  };
  class OptimiseArrangmentNode:public Node {
    float data_multiplier = 50.0;
    float smoothness_multiplier = 1.0;
    bool preset_labels = false;
    int n_iterations = 3;
    int graph_cut_impl = 2;

    float z_percentile = 0.9;

    public:
    using Node::Node;
    void init() override {
      add_input("arrangement", typeid(Arrangement_2));
      add_input("pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
      add_output("arrangement", typeid(Arrangement_2));

      add_param(ParamBoundedInt(graph_cut_impl, 0, 2, "graph_cut_impl", "Graph cut implementation"));
      add_param(ParamBoundedInt(n_iterations, 0, 100, "n_iterations", "Number of iterations"));
      add_param(ParamBoundedFloat(data_multiplier, 0.001, 100, "data_multiplier", "Multiplier on data term"));
      add_param(ParamBoundedFloat(smoothness_multiplier, 0.001, 100, "smoothness_multiplier", "Multiplier on smoothness term"));
      add_param(ParamBool(preset_labels, "preset_labels", "Preset face labels"));
      add_param(ParamBoundedFloat(z_percentile, 0, 1, "z_percentile",  "Elevation percentile"));
    }
    void process() override;
  };

  class ArrDissolveNode: public Node {
    bool dissolve_seg_edges = true;
    bool dissolve_step_edges = false;
    bool dissolve_outside_fp = true;
    bool dissolve_all_interior = false;
    bool skip_execution = false;
    // bool remove_duplicates = true;
    float step_height_threshold = 1.0;
    // int dupe_threshold_exp = 3;
    public:
    using Node::Node;
    void init() override {
      add_input("arrangement", typeid(Arrangement_2));
      add_input("heightfield", typeid(RasterTools::Raster));
      
      add_output("global_elevation_97p", typeid(float));
      add_output("global_elevation_70p", typeid(float));
      add_output("global_elevation_50p", typeid(float));
      add_output("global_elevation_min", typeid(float));
      add_output("global_elevation_max", typeid(float));

      add_output("arrangement", typeid(Arrangement_2));

      add_param(ParamBool(dissolve_outside_fp, "dissolve_outside_fp", "Dissolve edges outside footprint"));
      add_param(ParamBool(dissolve_all_interior, "dissolve_all_interior", "Dissolve all edges in intrerior of footprint"));
      add_param(ParamBool(dissolve_seg_edges, "dissolve_seg_edges", "Dissolve same label cells"));
      add_param(ParamBool(dissolve_step_edges, "dissolve_step_edges", "Dissolve step edges (check for LoD1.3)"));
      add_param(ParamBool(skip_execution, "skip_execution", "Switch of autorun prior to execution in run_all()"));
      // add_param(ParamBool(remove_duplicates, "remove_duplicates", "Remove duplicates"));
      add_param(ParamBoundedFloat(step_height_threshold, 0, 10, "step_height_threshold",  "step_height_threshold"));
      // add_param(ParamBoundedInt(dupe_threshold_exp, 0, 15, "dupe_threshold_exp", "10-base exponent to set duplication threshold. Eg a value of 2 yields a value of 10^(-2) = 0.01"));
    }
    void process() override;
    bool parameters_valid() override {
      return !skip_execution;
    }
  };

  class LinearRingtoRingsNode:public Node {
    public:
    using Node::Node;
    void init() override {
      add_input("linear_ring", typeid(LinearRing));
      add_output("linear_rings", typeid(LinearRingCollection));
    }
    void process() override;
  };

  class DetectLinesNode:public Node {
    float dist_thres = 0.4;
    std::pair<int,int> min_cnt_range = {5,10};
    int min_cnt_range_lower = 5;
    int min_cnt_range_upper = 10;
    int k = 10;
    float snap_threshold = 1;
    float line_extend = 0.05;
    bool perform_chaining = true;
    bool remove_overlap = true;
    public:
    typedef std::pair<size_t,size_t> IDPair;
    struct Cmp {
      bool operator()(const IDPair& lhs, const IDPair& rhs) const { 
          return lhs.first < rhs.first; 
      }
    };
    typedef std::map<IDPair, size_t, Cmp> RingSegMap;

    using Node::Node;
    void init() override {
      add_vector_input("edge_points", {typeid(LinearRing)});
      add_input("roofplane_ids", typeid(vec1i));
      add_input("pts_per_roofplane", typeid(IndexedPlanesWithPoints ));
      add_output("edge_segments", typeid(SegmentCollection));
      add_output("lines3d", typeid(SegmentCollection));
      add_output("ring_edges", typeid(SegmentCollection));
      add_output("ring_idx", typeid(std::unordered_map<size_t,std::vector<size_t>>));
      add_output("ring_id", typeid(vec1i));
      add_output("ring_order", typeid(vec1i));
      add_output("is_start", typeid(vec1i));

      add_param(ParamFloat(dist_thres, "dist_thres", "dist_thres"));
      add_param(ParamInt(min_cnt_range_lower, "min_cnt_range_lower", "Minimum count lower"));
      add_param(ParamInt(min_cnt_range_upper, "min_cnt_range_upper", "Minimum count upper"));
      add_param(ParamInt(k, "k", "k"));
      add_param(ParamFloat(snap_threshold, "snap_threshold", "Chain snap thres"));
      add_param(ParamFloat(line_extend, "line_extend", "Extend lines"));
      add_param(ParamBool(perform_chaining, "perform_chaining", "Perform chaining"));
      add_param(ParamBool(remove_overlap, "remove_overlap", "Remove overlap"));
    }
    inline void detect_lines_ring_m1(linedect::LineDetector& LD, SegmentCollection& segments_out);
    inline size_t detect_lines_ring_m2(linedect::LineDetector& LD, Plane&, SegmentCollection& segments_out);
    inline void detect_lines(linedect::LineDetector& LD);
    void process() override;
  };

  // class ClassifyEdgePointsNode:public Node {
  //   public:
  //   using Node::Node;
  //   void init() override {
  //     add_output("edge_points", typeid(std::vector<linedect::Point>));
  //     add_output("edge_points_vec3f", typeid(vec3f));
  //     add_input("points", typeid(PNL_vector));

  //     add_param("classify_jump_count_min", (int) 1);
  //     add_param("classify_jump_count_max", (int) 5);
  //     add_param("classify_line_dist", (float) 0.005);
  //     add_param("classify_jump_ele", (float) 1.0);
  //   }

  //   void gui(){
  //     ImGui::InputInt("Jump cnt min", &param<int>("classify_jump_count_min"));
  //     ImGui::InputInt("Jump cnt max", &param<int>("classify_jump_count_max"));
  //     ImGui::InputFloat("Line dist", &param<float>("classify_line_dist"), 0.01, 1);
  //     ImGui::InputFloat("Elevation jump", &param<float>("classify_jump_ele"), 0.01, 1);
  //   }

  //   void process() override;
  // };

  class DetectPlanesNode:public Node {
    float horiz_min_count = 0.95;
    int metrics_normal_k = 5;
    int metrics_plane_k = 15;
    int metrics_plane_min_points = 20;
    float metrics_plane_epsilon = 0.2;
    float metrics_plane_normal_threshold = 0.75;
    float metrics_is_horizontal_threshold = 0.97;
    float metrics_probability_ransac = 0.05;
    float metrics_cluster_epsilon_ransac = 0.3;
    float metrics_is_wall_threshold = 0.3;
    int n_refit = 5;
    bool use_ransac = false;
    // float roof_percentile=0.5;

    // plane regularisation
    float maximum_angle_ = 25;
    float maximum_offset_ = 0.5;
    bool regularize_parallelism_ = false;
    bool regularize_orthogonality_ = false;
    bool regularize_coplanarity_ = false;
    bool regularize_axis_symmetry_ = false;

    public:
    using Node::Node;
    void init() override {
      add_input("points", typeid(PointCollection));
      add_output("plane_id", typeid(vec1i));
      add_output("is_wall", typeid(vec1i));
      add_output("is_horizontal", typeid(vec1i));
      
      add_output("pts_per_roofplane", typeid(IndexedPlanesWithPoints ));

      add_output("roof_pt_cnt", typeid(int));
      add_output("wall_pt_cnt", typeid(int));
      add_output("unsegmented_pt_cnt", typeid(int));
      add_output("roof_type", typeid(std::string));
      add_output("roof_elevation_70p", typeid(float));
      add_output("roof_elevation_50p", typeid(float));
      add_output("roof_elevation_min", typeid(float));
      add_output("roof_elevation_max", typeid(float));
      add_output("horiz_roofplane_cnt", typeid(float));
      add_output("slant_roofplane_cnt", typeid(float));
      add_output("total_roofplane_cnt", typeid(int));
      add_output("plane_adj", typeid(std::map<size_t, std::map<size_t, size_t>>));

      add_param(ParamBool(use_ransac, "use_ransac", "Use ransac instead of region growing plane detection"));
      add_param(ParamFloat(horiz_min_count, "horiz_min_count", "Mininmal point count for horizontal planes"));
      add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
      add_param(ParamInt(metrics_plane_k, "metrics_plane_k", "Number of neighbours used during region growing plane detection"));
      add_param(ParamInt(metrics_plane_min_points, "metrics_plane_min_points", "Minimum number of points in a plane"));
      add_param(ParamFloat(metrics_plane_epsilon, "metrics_plane_epsilon", "Plane epsilon"));
      add_param(ParamFloat(metrics_cluster_epsilon_ransac, "metrics_cluster_epsilon_ransac", "Cluster epsilon RANSAC only"));
      add_param(ParamFloat(metrics_probability_ransac, "metrics_probability_ransac", "Probability RANSAC only"));
      add_param(ParamFloat(metrics_plane_normal_threshold, "metrics_plane_normal_threshold", "Plane normal angle threshold"));
      add_param(ParamFloat(metrics_is_horizontal_threshold, "metrics_is_horizontal_threshold", "Threshold for horizontal plane detection (expressed as angle wrt unit verctor in +z direction)"));
      add_param(ParamFloat(metrics_is_wall_threshold, "metrics_is_wall_threshold", "Wall angle thres"));
      add_param(ParamInt(n_refit, "n_refit", "Refit every n points"));

      add_param(ParamBool(regularize_parallelism_, "regularize_parallelism", "regularizeparallelism"));
      add_param(ParamBool(regularize_orthogonality_, "regularize_orthogonality", "regularizeorthogonality"));
      add_param(ParamBool(regularize_coplanarity_, "regularize_coplanarity", "regularizecoplanarity"));
      add_param(ParamBool(regularize_axis_symmetry_, "regularize_axis_symmetry", "regularize_axis_symmetry"));
      add_param(ParamFloat(maximum_angle_, "maximum_angle", "Plane regularisation: maximum allowed angle in degrees between plane normals used for parallelism, orthogonality, and axis symmetry"));
      add_param(ParamFloat(maximum_offset_, "maximum_offset", "Plane regularisation: maximum allowed orthogonal distance between two parallel planes such that they are considered to be coplanar"));
      // add_param(ParamBoundedFloat(roof_percentile, 0, 1, "roof_percentile",  "Roof elevation percentile"));
    }
    // void before_gui(){
    //   auto param_count = std::get<ParamFloat>(parameters.at("horiz_min_count"));
    //   auto param_ishoriz = std::get<ParamFloat>(parameters.at("metrics_is_horizontal_threshold"));
    //   param_count.set_visible(only_horizontal);
    //   param_ishoriz.set_visible(only_horizontal);
    // }
    void process() override;
  };

  class LASInPolygonsNode:public Node {
    std::string filepaths_ = "";
    float cellsize = 50.0;
    float buffer = 1.0;
    float ground_percentile=0.05;
    float max_density_delta=0.05;
    float coverage_threshold=2.0;
    int ground_class = 2;
    int building_class = 6;
    bool clear_if_insufficient = true;
    std::string wkt_="";
    public:
    using Node::Node;
    void init() override {
      add_vector_input("polygons", typeid(LinearRing));
      add_vector_input("buf_polygons", typeid(LinearRing));
      add_vector_output("point_clouds", typeid(PointCollection));
      add_vector_output("ground_elevations", typeid(float));
      add_vector_output("poly_areas", typeid(float));
      add_vector_output("poly_pt_counts_bld", typeid(int));
      add_vector_output("poly_pt_counts_grd", typeid(int));
      add_vector_output("poly_ptcoverage_class", typeid(std::string));
      add_vector_output("poly_densities", typeid(float));
      add_param(ParamInt(ground_class, "ground_class", "LAS class number to use for ground"));
      add_param(ParamInt(building_class, "building_class", "LAS class number to use for buildings"));

      add_param(ParamPath(filepaths_, "las_filepaths", "LAS filepaths"));
      add_param(ParamBoundedFloat(cellsize, 1, 1000, "cellsize",  "Grid index cellsize"));
      add_param(ParamBoundedFloat(buffer, 0.1, 100, "buffer", "Query buffer"));
      add_param(ParamBoundedFloat(ground_percentile, 0, 1, "ground_percentile",  "Ground elevation percentile"));
      add_param(ParamBoundedFloat(max_density_delta, 0, 1, "Max_density_delta",  "Used for deciding to what footprint to assign a point that is inside multiple footprints. If the difference in point densities is higher than this threshold we pick the candidate footprint with the highest point density. Otherwise we pick the footprint with the highest average elevation."));
      add_param(ParamBoundedFloat(coverage_threshold, 0.5, 3, "coverage_threshold",  "Threshold to classify if a footprint has sufficient building point_coverage. If the building point density is below this threshold times the of standard deviation below the mean point density, the coverage is classified as insufficient."));
      add_param(ParamBool(clear_if_insufficient, "clear_if_insufficient",  "Clear point clouds for footprints with insufficient coverage"));
      add_param(ParamString(wkt_, "wkt",  "Override CRS"));
    }
    void process() override;
  };

#ifdef GFP_WITH_PDAL
  class EptInPolygonsNode:public Node {
    std::string dirpath = "";
    std::string filter_limits = "Classification[2:2],Classification[6:6]";
    float cellsize = 50.0;
    float buffer = 1.0;
    float ground_percentile=0.05;
  public:
    using Node::Node;
    void init() override {
      add_vector_input("polygons", typeid(LinearRing));
      add_vector_input("buf_polygons", typeid(LinearRing));
      add_vector_output("point_clouds", typeid(PointCollection));
      add_vector_output("ground_point_clouds", typeid(PointCollection));

      add_param(ParamPath(dirpath, "dirpath", "EPT directory"));
      // add_param("filter_limits", ParamString(filter_limits, "PDAL Range filter"));
      add_param(ParamBoundedFloat(cellsize, 1, 1000, "cellsize",  "Grid index cellsize"));
      add_param(ParamBoundedFloat(buffer, 0.1, 100, "buffer", "Query buffer"));
      add_param(ParamBoundedFloat(ground_percentile, 0, 1, "ground_percentile",  "Ground elevation percentile"));
    }
    void process() override;
  };
#endif
  
  class BuildingSelectorNode:public Node {
    public:
    int building_id=0, polygon_count=0;
    using Node::Node;
    void init() override {
      add_vector_input("point_clouds", typeid(PointCollection));
      add_vector_input("ground_point_clouds", typeid(PointCollection));
      add_vector_input("polygons", typeid(LinearRing));
      add_vector_input("ground_elevations", typeid(float));
      add_output("point_cloud", typeid(PointCollection));
      add_output("ground_point_cloud", typeid(PointCollection));
      add_output("polygon", typeid(LinearRing));
      add_output("ground_elevation", typeid(float));

      add_param(ParamBoundedInt(building_id, 0, polygon_count-1,  "building_id", "building_id"));
    }
    void on_receive(gfInputTerminal& it) {
      if (it.get_name() == "polygons") {
        polygon_count = vector_input("polygons").size();
        auto param = static_cast<ParamBoundedInt*>(parameters.at("building_id").get());
        param->set_bounds(0, polygon_count-1);
      }
    }
    void process() override;
  };

  class RegulariseLinesNode:public Node {
    float dist_threshold = 0.5;
    float angle_threshold = 0.15;
    float extension = 1.0;
    bool merge_intersection_lines = false;

    public:
    using Node::Node;
    void init() override {
      add_input("edge_segments", {typeid(SegmentCollection), typeid(LineString)});
      add_input("ints_segments", typeid(SegmentCollection));
      // add_input("footprint", typeid(LinearRing));

      add_vector_output("regularised_edges", typeid(Segment));
      add_vector_output("exact_regularised_edges", typeid(linereg::Segment_2));
      add_output("edges_out_", typeid(SegmentCollection));
      add_output("priorities", typeid(vec1i));
      add_output("angle_cluster_id", typeid(vec1i));
      add_output("dist_cluster_id", typeid(vec1i));
      add_output("exact_footprint_out", typeid(linereg::Polygon_with_holes_2));
      add_output("n_angle_clusters", typeid(int));
      
      add_param(ParamFloat(dist_threshold, "dist_threshold", "Distance threshold"));
      add_param(ParamFloat(extension, "extension", "Line extension after regularisation"));
      add_param(ParamBoundedFloat(angle_threshold, 0.01, 3.1415, "angle_threshold", "Angle threshold"));
      add_param(ParamBool(merge_intersection_lines, "merge_intersection_lines", "merge_intersection_lines"));

    }
    void process() override;
  };

  class ClusterLinesNode:public Node {
    float dist_threshold = 0.5;
    float angle_threshold = 0.15;

    public:
    using Node::Node;
    void init() override {
      add_input("segments", typeid(SegmentCollection));
      add_output("segments", typeid(SegmentCollection));
     
      add_param(ParamFloat(dist_threshold, "dist_threshold", "Distance threshold"));
      add_param(ParamBoundedFloat(angle_threshold, 0.01, 3.1415, "angle_threshold", "Angle threshold"));

    }
    void process() override;
  };

  class RegulariseRingsNode:public Node {
    float dist_threshold = 0.5;
    float angle_threshold = 0.15;
    float snap_threshold = 1.0;
    bool weighted_avg = false;
    bool angle_per_distcluster = false;
    bool regularise_fp = false;
    bool recompute_rings = false;
    float fp_offset = 0.01;

    public:
    using Node::Node;
    void init() override {
      add_input("edge_segments", typeid(SegmentCollection));
      add_input("ints_segments", typeid(SegmentCollection));
      add_input("ring_idx", typeid(std::unordered_map<size_t,std::vector<size_t>>));
      // add_input("ring_id", typeid(vec1i));
      // add_input("ring_order", typeid(vec1i));
      // add_input("edge_segments", typeid(SegmentCollection));
      add_input("footprint", typeid(LinearRing));
      add_vector_output("edges_out", typeid(Segment));
      add_output("edges_out_", typeid(SegmentCollection));
      add_output("priorities", typeid(vec1i));
      add_output("angle_cluster_id", typeid(vec1i));
      add_output("dist_cluster_id", typeid(vec1i));
      // add_output("rings_out", typeid(LinearRingCollection));
      // add_output("footprint_out", typeid(LinearRing));
      add_output("rings_out", typeid(LinearRingCollection));
      add_output("plane_id", typeid(vec1i));
      add_output("exact_rings_out", typeid(std::unordered_map<size_t, linereg::Polygon_2>));
      add_output("exact_footprint_out", typeid(linereg::Polygon_with_holes_2));
      // add_output("footprint_labels", typeid(vec1i));
      // add_output("line_clusters", TT_any); // ie a LineCluster
      // add_output("tmp_vec3f", typeid(vec3f));
      add_param(ParamFloat(dist_threshold, "dist_threshold", "Distance threshold"));
      add_param(ParamBoundedFloat(angle_threshold, 0.01, 3.1415, "angle_threshold", "Angle threshold"));
      add_param(ParamBoundedFloat(snap_threshold, 0.01, 10, "snap_threshold", "Snap threshold"));
      add_param(ParamBool(weighted_avg, "weighted_avg", "weighted_avg"));
      add_param(ParamBool(angle_per_distcluster, "angle_per_distcluster", "angle_per_distcluster"));
      add_param(ParamBool(regularise_fp, "regularise_fp", "regularise_fp"));
      add_param(ParamBool(recompute_rings, "recompute_rings", "recompute_rings"));
      add_param(ParamBoundedFloat(fp_offset, 0.01, 10, "fp_offset", "fp_offset"));
    }
    void process() override;
  };


  class PrintResultNode:public Node {
    public:
    using Node::Node;
    void init() override {
      add_input("in", {typeid(float)});
    }
    std::string info() override {
      return "Result: " + std::to_string(input("in").get<float>());
    }
    void process() override {};
  };


  class PlaneIntersectNode:public Node {
    int min_neighb_pts = 5;
    float min_dist_to_line = 1.0;
    float min_length = 0;
    public:
    using Node::Node;
    void init() override {
      add_input("pts_per_roofplane", 
        {typeid(IndexedPlanesWithPoints)});
      add_input("plane_adj", 
        {typeid(std::map<size_t, std::map<size_t, size_t>>)});
      // add_vector_input("alpha_rings", 
        // {typeid(LinearRing)});

      add_output("lines", typeid(SegmentCollection));

      add_param(ParamInt(min_neighb_pts, "min_neighb_pts", "Minimum number of neighbouring points"));
      add_param(ParamFloat(min_dist_to_line, "min_dist_to_line", "Minimum number of neighbouring points"));
      add_param(ParamFloat(min_length, "min_length", "Minimum length of segment in order to be outputted"));
    }
    void process() override;
  };


  class SimplifyPolygonNode:public Node {
    float threshold_stop_cost = 0.005;
    public:
    using Node::Node;
    void init() override {
      add_vector_input("polygons", typeid(LinearRing));
      add_vector_output("polygons_simp", typeid(LinearRing));
      // add_output("polygon_simp", typeid(LinearRing));

      add_param(ParamBoundedFloat(threshold_stop_cost, 0, 1000, "threshold_stop_cost",  "threshold_stop_cost"));
    }
    // void on_change_parameter(std::string name, ParameterVariant& param){
    //   if (name == "threshold_stop_cost") {
    //     manager.run(*this);
    //   }
    // }
    void process() override;
  };


  class PC2MeshQualityNode:public Node {

    public:
    using Node::Node;
    void init() override {
      add_input("ipoints", {typeid(PointCollection), typeid(IndexedPlanesWithPoints)});
      add_input("triangles", typeid(MultiTriangleCollection));
      add_input("face_ids", typeid(vec1i));
      
      add_output("point_errors", typeid(vec1f));
      add_output("face_errors", typeid(vec1f));
      add_output("mesh_error_f", typeid(float));
      add_output("error_hist", typeid(std::string));
      add_output("mesh_error", typeid(vec1f));
      add_output("m2pc_error_hist", typeid(std::string));
      add_output("m2pc_error_max", typeid(float));

    }
    void process() override;
  };

  class PCRasteriseNode:public Node {
    float cellsize = 0.5;
    bool use_footprint_ = false;

    public:
    using Node::Node;
    void init() override {
      add_input("points", typeid(PointCollection));
      add_input("footprint", typeid(LinearRing));
      
      add_output("heightfield", typeid(RasterTools::Raster));
      add_output("grid_points", typeid(PointCollection));
      add_output("values", typeid(vec1f));
      add_output("image_fp", typeid(geoflow::Image));
      add_output("image_max", typeid(geoflow::Image));
      add_output("image_min", typeid(geoflow::Image));
      add_output("image_cnt", typeid(geoflow::Image));
      add_output("image_med", typeid(geoflow::Image));
      add_output("image_avg", typeid(geoflow::Image));
      add_output("image_var", typeid(geoflow::Image));

      add_param(ParamBoundedFloat(cellsize, 0, 50, "cellsize",  "cellsize"));
      add_param(ParamBool(use_footprint_, "use_footprint",  "use_footprint for setting raster bounds."));
    }
    void process() override;

    bool inputs_valid() override {
      for (auto& [name,iT] : input_terminals) {
        if (name == "footprint" && !use_footprint_) {
          continue;
        } else if (!iT->has_data()) {
          return false;
        }
      }
      return true;
    }
  };

  class RoofPartitionRasteriseNode:public Node {
    float cellsize = 0.5;
    bool use_planes_ = true;

    public:
    using Node::Node;
    void init() override {
      add_input("arrangement", typeid(Arrangement_2));
      add_input("footprint", typeid(LinearRing));
      add_input("h_ground", typeid(float));

      // add_output("values", typeid(vec1f));
      add_output("image", typeid(geoflow::Image));
      add_output("volume", typeid(float));

      add_param(ParamBool(use_planes_, "use_planes",  "Use actual planes to determine elevation of pixels instead of 70p elevation of whole roofpart."));
      add_param(ParamBoundedFloat(cellsize, 0, 50, "cellsize",  "cellsize"));
    }
    void process() override;
    void rasterise_arrangement(Arrangement_2& arr, RasterTools::Raster& r, size_t& data_pixel_cnt, bool use_planes);

  };

  class RoofPartition3DBAGRasteriseNode:public Node {
    float cellsize = 0.5;

    public:
    using Node::Node;
    void init() override {
      add_input("lod12_roofparts", typeid(geoflow::LinearRing));
      add_input("lod13_roofparts", typeid(geoflow::LinearRing));
      add_input("lod22_roofparts", typeid(geoflow::LinearRing));

      // add_output("values", typeid(vec1f));
      add_poly_output("lod12_hattr", {typeid(float)});
      add_poly_output("lod13_hattr", {typeid(float)});
      add_poly_output("lod22_hattr", {typeid(float)});

      add_param(ParamBoundedFloat(cellsize, 0, 50, "cellsize",  "cellsize"));
    }
    void process() override;
  };

  class BuildingRasteriseNode:public Node {
    float cellsize = 0.5;
    bool use_tin = false;
    bool use_ground_pts = true;

    float angle_thres = 5;
    float normal_angle_thres = 10;
    float area_thres = 0.3;
    float len_thres = 0.65;
    bool do_angle_thres = false;
    bool do_normal_angle_thres = false;
    bool do_area_thres = false;
    bool do_len_thres = true;

    public:
    using Node::Node;
    void init() override {
      add_input("points", typeid(PointCollection));
      add_input("ground_points", typeid(PointCollection));
      add_input("h_ground", typeid(float));
      add_input("footprint", typeid(LinearRing));
      
      add_output("heightfield", typeid(RasterTools::Raster));
      add_output("grid_points", typeid(PointCollection));
      add_output("values", typeid(vec1f));
      add_output("image", typeid(geoflow::Image));

      add_param(ParamBoundedFloat(cellsize, 0, 50, "cellsize",  "Cellsize"));
      add_param(ParamBool(use_tin, "use_tin",  "Use TIN for rasterisation"));
      add_param(ParamBool(use_ground_pts, "use_ground_pts",  "Rasterise the ground points"));

      add_param(ParamBoundedFloat(angle_thres, 0, 180, "angle_thres",  "angle_thres"));
      add_param(ParamBoundedFloat(normal_angle_thres, 0, 180, "normal_angle_thres",  "normal_angle_thres"));
      add_param(ParamBoundedFloat(len_thres, 0, 3, "len_thres",  "len_thres"));
      add_param(ParamBoundedFloat(area_thres, 0, 100, "area_thres",  "area_thres"));
      add_param(ParamBool(do_angle_thres, "do_angle_thres",  "do_angle_thres"));
      add_param(ParamBool(do_normal_angle_thres, "do_normal_angle_thres",  "do_normal_angle_thres"));
      add_param(ParamBool(do_len_thres, "do_len_thres",  "do_len_thres"));
      add_param(ParamBool(do_area_thres, "do_area_thres",  "do_area_thres"));
    }
    void process() override;
  };

  class SegmentRasteriseNode:public Node {
    float cellsize = 0.05;
    float thres_alpha = 0.25;
    bool use_ground = false;
    int megapixel_limit = 600;
    bool fill_nodata_ = true;
    int fill_nodata_window_size_ = 5;
    public:
    using Node::Node;
    void init() override {
      // add_vector_input("alpha_rings", typeid(LinearRing));
      add_vector_input("triangles", typeid(TriangleCollection));
      add_vector_input("ground_triangles", typeid(TriangleCollection));
      add_vector_input("alpha_rings", typeid(LinearRing));
      add_input("roofplane_ids", typeid(vec1i));
      add_input("pts_per_roofplane", typeid(IndexedPlanesWithPoints));
      // add_input("heightfield", typeid(RasterTools::Raster));
      add_output("heightfield", typeid(RasterTools::Raster));
      
      add_output("heightfield", typeid(RasterTools::Raster));
      add_output("grid_points", typeid(PointCollection));
      add_output("values", typeid(vec1f));
      add_output("data_area", typeid(float));

      add_param(ParamBoundedFloat(cellsize, 0, 50, "cellsize",  "cellsize"));
      add_param(ParamBool(use_ground, "use_ground",  "Rasterise the ground_triangles input"));
      add_param(ParamBoundedInt(fill_nodata_window_size_, 0, 50, "fill_nodata_window_size",  "fill_nodata_window_size"));
      add_param(ParamBool(fill_nodata_, "fill_nodata",  "Fill nodata values with NN interpolation"));
      add_param(ParamInt(megapixel_limit, "megapixel_limit", "Max size of raster in megalpixels. If exceeded cellsize is increased."));

    }
    bool inputs_valid() override {
      for (auto& [name,iT] : input_terminals) {
        if (name == "ground_triangles" && !use_ground) {
          continue;
        } else if (!iT->is_touched()) {
          return false;
        }
      }
      return true;
    }
    void rasterise_input(gfSingleFeatureInputTerminal& input_triangles, RasterTools::Raster& r, size_t& data_pixel_cnt);
    void process() override;
  };

  // class GridMaxNode:public Node {
  //   public:
  //   using Node::Node;
  //   void init() override {
  //     add_input("grid_1", typeid(RasterTools::Raster));
  //     add_input("grid_2", typeid(RasterTools::Raster));
      
  //     add_output("heightfield", typeid(RasterTools::Raster));
  //     add_output("grid_points", typeid(PointCollection));
  //     // add_output("values", typeid(vec1f));
  //   }
  //   void process() override;
  // };

  class PolygonUnionNode:public Node {
    public:
    using Node::Node;
    void init() override {
      add_vector_input("polygons", typeid(LinearRing));
      add_vector_output("polygons", typeid(LinearRing));
      add_vector_output("holes", typeid(LinearRing));
    }
    void process() override;
  };

  class Filter25DNode:public Node {
    float cellsize = 0.1;
    float angle_thres = 5;
    float normal_angle_thres = 10;
    float area_thres = 0.3;
    float len_thres = 0.65;
    bool do_angle_thres = false;
    bool do_normal_angle_thres = false;
    bool do_area_thres = false;
    bool do_len_thres = true;
    // int count_thres = 1;
  public:
    using Node::Node;
    template<typename T> void mark_big_triangles(T& dt);
    void init() override {
      add_input("points", {typeid(PointCollection), typeid(IndexedPlanesWithPoints)});
      add_output("points", typeid(PointCollection));
      add_output("points_filtered", typeid(PointCollection));
      add_output("triangles", typeid(TriangleCollection));
      add_output("heightfield", typeid(RasterTools::Raster));
      add_output("heightfield_pts", typeid(PointCollection));

      add_param(ParamBoundedFloat(cellsize, 0, 50, "cellsize",  "cellsize"));
      // add_param(ParamBoundedInt(count_thres, 0, 10, "count_thres", "count_thres"));
      add_param(ParamBoundedFloat(angle_thres, 0, 180, "angle_thres",  "angle_thres"));
      add_param(ParamBoundedFloat(normal_angle_thres, 0, 180, "normal_angle_thres",  "normal_angle_thres"));
      add_param(ParamBoundedFloat(len_thres, 0, 3, "len_thres",  "len_thres"));
      add_param(ParamBoundedFloat(area_thres, 0, 100, "area_thres",  "area_thres"));
      add_param(ParamBool(do_angle_thres, "do_angle_thres",  "do_angle_thres"));
      add_param(ParamBool(do_normal_angle_thres, "do_normal_angle_thres",  "do_normal_angle_thres"));
      add_param(ParamBool(do_len_thres, "do_len_thres",  "do_len_thres"));
      add_param(ParamBool(do_area_thres, "do_area_thres",  "do_area_thres"));
    }
    void process() override;
  };

  class PolygonTriangulatorNode:public Node {
    int dupe_threshold_exp = 6;
    bool output_all_triangles = false;
    bool output_mtc_for_every_input = false;

    void triangulate_polygon(LinearRing& ring, vec3f& normals, TriangleCollection& triangles, size_t& ring_id, vec1i& ring_ids);
    public:
    using Node::Node;
    void init() override {
      add_vector_input("polygons", {typeid(LinearRing), typeid(std::unordered_map<int, Mesh>), typeid(Mesh)});
      // add_vector_output("dupe_rings", typeid(LinearRing));
      add_output("triangles", typeid(TriangleCollection));
      add_output("multi_triangle_collections", typeid(MultiTriangleCollection));
      add_output("normals", typeid(vec3f));
      add_output("volumes", typeid(float));
      add_output("ring_ids", typeid(vec1i));
      // add_output("nesting_levels", typeid(vec1i));

      // add_output("edges", typeid(SegmentCollection));
      // add_output("edges_constr", typeid(vec1i));

      add_param(ParamBool(output_all_triangles, "output_all_triangles",  "Also output triangles in holes and outside convex hull."));
      add_param(ParamBool(output_mtc_for_every_input, "output_mtc_for_every_input",  "Output a MultiTriangleCollection for every input. Otherwise aggregate all inputs into one MultiTriangleCollection. Only applies to Mesh type inputs."));
      add_param(ParamInt(dupe_threshold_exp, "dupe_threshold_exp", "Dupe tolerance exponent"));

    }
    void process() override;
  };

  // convert from dreadful multitrianglecollection to slightly less dreadful meshmap.
  class MTC2MMNode:public Node {
    public:
    using Node::Node;
    void init() override {
      add_input("multi_triangle_collections", typeid(MultiTriangleCollection));
      add_output("meshmap", typeid(std::unordered_map<int, Mesh>));

    }
    void process() override;
  };

  class FacesSelectorNode:public Node {
    public:
    using Node::Node;
    bool inputs_valid() override {
      auto& skip_term = input("skip");
      auto& replace_term = input("replace");
      if ((!skip_term.has_data() || (!replace_term.has_data()))) {
        return false;
      } else {
        bool skip = skip_term.get<bool>();
        if (skip) {
          return vector_input("faces_A").has_data();
        } else {
          return vector_input("faces_B").has_data() && vector_input("faces_A").has_data();
        }
      }
    }
    void init() override {
      add_input("skip", typeid(bool));
      add_input("replace", typeid(bool));
      add_vector_input("faces_A", typeid(std::unordered_map<int, Mesh>));
      add_vector_input("faces_B", typeid(std::unordered_map<int, Mesh>));
      add_vector_output("faces", typeid(std::unordered_map<int, Mesh>));
    }
    void process() override {
      auto skip = input("skip").get<bool>();
      auto replace = input("replace").get<bool>();
      if (skip) {
        if (replace)
          vector_output("faces") = vector_input("faces_A").get_data_vec();
      } else {
        vector_output("faces") = vector_input("faces_B").get_data_vec();
      }
    };
  };


  class AttributeTesterNode:public Node {
    std::string attribute_name;
    public:
    using Node::Node;
    void init() override {
      add_poly_input("attributes", {typeid(bool), typeid(int), typeid(float), typeid(std::string), typeid(std::string), typeid(Date), typeid(Time), typeid(DateTime)});
      add_output("result", typeid(bool));

      add_param(ParamString(attribute_name, "attribute_name", "Attribute name (should be a boolean attribute)"));

    }
    void process() override {
      bool result = false;
      for (auto &iterm : poly_input("attributes").sub_terminals()) {
        if (iterm->get_name() == manager.substitute_globals(attribute_name)) {
          if (iterm->accepts_type(typeid(bool))) {
            // std::cout << "Detected attribute of type bool\n";
            result = iterm->get<bool>();
          } else {
            // std::cout << "Attribute type is not bool\n";
          }
        }
      }
      output("result").set(result);
    };
  };

  class SkipReplaceTesterNode:public Node {
    std::string attribute_name;
    public:
    using Node::Node;
    void init() override {
      add_input("skip", typeid(bool));
      add_vector_input("alpha_rings", {typeid(LinearRing)});
      add_input("roof_type", typeid(std::string));
      add_output("roof_type", typeid(std::string));
      add_output("skip", typeid(bool));
      add_output("replace", typeid(bool));

      add_param(ParamString(attribute_name, "attribute_name", "Attribute name (should be a boolean attribute)"));

    }
    bool inputs_valid() override {
      auto& skip_term = input("skip");
      auto& ar_term = vector_input("alpha_rings");
      if(input("roof_type").has_data() && skip_term.has_data()) {
        auto skip = skip_term.get<bool>();
        if (skip) return true;
        else return ar_term.is_touched();
      } 
      return false;
      
    }
    void process() override {
      bool skip = input("skip").get<bool>();
      bool replace = skip;
      
      auto roof_type = input("roof_type").get<std::string>();
      if(vector_input("alpha_rings").is_touched() && vector_input("alpha_rings").size()==0 && roof_type != "no points" && roof_type != "no planes") {
        roof_type = "no points";
      }
      if (roof_type == "no points" || roof_type == "no planes") {{
        skip=true;
        replace=false;
      }}

      output("skip").set(skip);
      output("replace").set(replace);
      output("roof_type").set(roof_type);
    };
  };

  class AttrRingsSelectorNode:public Node {
    public:
    using Node::Node;
    bool inputs_valid() override {
      auto& skip_term = input("skip");
      auto& replace_term = input("replace");
      if ((!skip_term.has_data() || (!replace_term.has_data()))) {
        return false;
      } else {
        bool skip = skip_term.get<bool>();
        bool replace = replace_term.get<bool>();
        if (skip) {
          if (replace)
            return vector_input("linear_rings_A").has_data();
          else 
            return true;
        } else {
          return vector_input("linear_rings_B").has_data() && vector_input("linear_rings_A").has_data();
        }
      }
    }
    void init() override {
      add_input("skip", typeid(bool));
      add_input("replace", typeid(bool));
      add_vector_input("linear_rings_A", typeid(LinearRing));
      add_vector_input("linear_rings_B", typeid(LinearRing));
      add_poly_input("attributes_A", {typeid(bool), typeid(int), typeid(float), typeid(std::string), typeid(std::string), typeid(Date), typeid(Time), typeid(DateTime)});
      add_poly_input("attributes_B", {typeid(bool), typeid(int), typeid(float), typeid(std::string), typeid(std::string), typeid(Date), typeid(Time), typeid(DateTime)});

      add_vector_output("linear_rings", typeid(LinearRing));
      add_poly_output("attributes", {typeid(bool), typeid(int), typeid(float), typeid(std::string), typeid(std::string), typeid(Date), typeid(Time), typeid(DateTime)});
    }
    void process() override {
      auto skip = input("skip").get<bool>();
      auto replace = input("replace").get<bool>();
      auto& out_attributes = poly_output("attributes");
      
      if (skip) {
        out_attributes = poly_input("attributes_A");

        auto &attr_datacov_term = out_attributes.add_vector("data_coverage", typeid(float));
        attr_datacov_term.push_back_any(std::any());
        auto &attr_ground_term = out_attributes.add_vector("is_ground", typeid(bool));
        attr_ground_term.push_back_any(std::any());
        auto &part_id_term = out_attributes.add_vector("building_part_id", typeid(int));
        part_id_term.push_back_any(std::any());
          
        if (replace) {
          vector_output("linear_rings") = vector_input("linear_rings_A").get_data_vec();
        } else{
          vector_output("linear_rings").push_back_any(std::any());
        }
        
      } else {
        vector_output("linear_rings") = vector_input("linear_rings_B").get_data_vec();
        out_attributes = poly_input("attributes_B");
      }
    };
  };

  class PolygonOffsetNode:public Node {
    public:
    using Node::Node;
    float offset = 4;

    void init() override {
      add_vector_input("polygons", typeid(LinearRing));
      add_vector_output("offset_polygons", typeid(LinearRing));

      add_param(ParamBoundedFloat(offset, -10, 10, "offset",  "offset"));   
    }
    void process() override;
  };

  class DataCoverageCalcNode:public Node {
    public:
    using Node::Node;

    void init() override {
      add_vector_input("ground_parts", typeid(LinearRing));
      add_input("footprint_polygon", typeid(LinearRing));
      add_input("data_area", typeid(float));
      
      add_output("data_coverage", typeid(float));
    }
    bool inputs_valid() override {
      return input("ground_parts").is_touched() && vector_input("footprint_polygon").has_data() && input("data_area").has_data();
    }
    void process() override;
  };

  class SnapRoundNode:public Node {
    float pixel_size = 0.001;

    public:
    using Node::Node;

    void init() override {
      add_input("arrangement", typeid(Arrangement_2));
      add_output("arrangement", typeid(Arrangement_2));

      add_param(ParamBoundedFloat(pixel_size, 0, 1, "pixel_size",  "pixel size"));
    }

    void process() override;
  };

  class TriSnapNode:public Node {
    // float pixel_size = 0.001;
    float dist_thres = 0.005;

    public:
    using Node::Node;

    void init() override {
      add_input("arrangement", typeid(Arrangement_2));
      add_output("arrangement", typeid(Arrangement_2));
      add_output("triangles_og", typeid(TriangleCollection));
      add_output("segment_ids_og", typeid(vec1i));
      add_output("triangles_snapped", typeid(TriangleCollection));
      add_output("segment_ids_snapped", typeid(vec1i));

      add_param(ParamBoundedFloat(dist_thres, 0, 1, "dist_thres",  "Snapping threshold"));
    }

    void process() override;
  };

  class SegmentExtendNode:public Node {
    // float pixel_size = 0.001;
    float extension = 0.005;

    public:
    using Node::Node;

    void init() override {
      add_input("segments", typeid(SegmentCollection));
      add_output("segments", typeid(SegmentCollection));

      add_param(ParamBoundedFloat(extension, 0, 50, "extension",  "Extension"));
    }

    void process() override;
  };

  class PointCloudMergerNode:public Node {

    public:
    using Node::Node;

    void init() override {
      add_poly_input("pointclouds", {typeid(PointCollection)});
      add_output("pointcloud", typeid(PointCollection));
    }

    void process() override;
  };

  class PC2PCDistancesCalculatorNode:public Node {

    public:
    using Node::Node;

    void init() override {
      add_input("pointcloud_a", {typeid(PointCollection)});
      add_input("pointcloud_b", {typeid(PointCollection)});
      add_output("errors_a_to_b", typeid(vec1f));
      add_output("errors_a_to_b_", typeid(float));
    }

    void process() override;
  };

  class PCFilterNode:public Node {
    
    float threshold = 1.0;

    public:
    using Node::Node;

    void init() override {
      add_input("pointcloud", {typeid(PointCollection)});
      add_input("values", {typeid(float)});
      add_output("pointcloud", typeid(PointCollection));

      add_param(ParamFloat(threshold, "threshold",  "Threshold "));
    }

    void process() override;
  };

  class PlaneIntersectAllNode:public Node {
    
    public:
    using Node::Node;

    void init() override {
      add_input("clipbbox", {typeid(PointCollection)});
      add_poly_input("planes", {typeid(Plane)});
      add_output("segments", typeid(SegmentCollection));
    }

    void process() override;
  };

  class RasterMergerNode:public Node {
    
    public:
    using Node::Node;

    void init() override {
      add_input("rasterbase", {typeid(RasterTools::Raster)});
      add_input("rasteradd", {typeid(RasterTools::Raster)});
      
      add_output("raster", typeid(RasterTools::Raster));
    }

    void process() override;
  };

  class ClusterPointCloudNode:public Node {
    bool flatten = true;
    float spacing = 1.0;
    // int metrics_normal_k = 5;

    public:
    using Node::Node;
    void init() override {
      add_input("points", typeid(PointCollection));
      add_output("cluster_id", typeid(vec1i));
      add_output("pts_per_roofplane", typeid(IndexedPlanesWithPoints));

      add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      add_param(ParamFloat(spacing, "spacing", "Radius threshold used in clustering algorithm"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

  class ContourRegulariserNode:public Node {
    // bool flatten = true;
    float min_length = 2.0;
    float max_angle = 20.0;
    float max_offset = 2.0;
    // int metrics_normal_k = 5;

    public:
    using Node::Node;
    void init() override {
      add_vector_input("polygons", typeid(LinearRing));
      // add_output("cluster_id", typeid(vec1i));
      add_vector_output("regularised_polygons", typeid(LinearRing));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      add_param(ParamFloat(min_length, "min_length", "Minimum length"));
      add_param(ParamFloat(max_angle, "max_angle", "Max angle"));
      add_param(ParamFloat(max_offset, "max_offset", "Max offset"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };


  class GridTilerNode : public Node {
    float cellsize_ = 500;

    public:
    using Node::Node;
    void init() override{
      add_vector_input("polygons", typeid(LinearRing));
      add_vector_output("tile_geoms", typeid(LinearRing));
      add_vector_output("tile_geom_ids", typeid(int));
      add_vector_output("tile_geom_cnts", typeid(int));
      add_vector_output("polygon_tile_ids", typeid(int));

      add_param(ParamFloat(cellsize_, "cellsize", "Cell size"));
    };
    void process() override {
      auto& polygons = input("polygons");
      auto& polygon_tile_ids = vector_output("polygon_tile_ids");
      auto& tile_geoms = vector_output("tile_geoms");
      auto& tile_geom_ids = vector_output("tile_geom_ids");
      auto& tile_geom_cnts = vector_output("tile_geom_cnts");

      Box bbox;
      for (size_t i=0; i<polygons.size(); ++i) {
        bbox.add(polygons.get<LinearRing&>(i).box());
      }
      auto grid = RasterTools::Raster(cellsize_, bbox.min()[0], bbox.max()[0], bbox.min()[1], bbox.max()[1]);

      std::unordered_map<int, int> tile_cnts;
      for (size_t i=0; i<polygons.size(); ++i) {
        auto c = polygons.get<LinearRing&>(i).box().center();
        int lc = (int) grid.getLinearCoord(c[0], c[1]);
        polygon_tile_ids.push_back(lc);
        tile_cnts[lc]++;
      }

      for(size_t col=0; col < grid.dimx_; ++col) {
        for(size_t row=0; row < grid.dimy_; ++row) {
          LinearRing g;
          g.push_back(arr3f{
            float(grid.minx_ + col*cellsize_),
            float(grid.miny_ + row*cellsize_),
            0            
          });
          g.push_back(arr3f{
            float(grid.minx_ + (col+1)*cellsize_),
            float(grid.miny_ + row*cellsize_),
            0            
          });
          g.push_back(arr3f{
            float(grid.minx_ + (col+1)*cellsize_),
            float(grid.miny_ + (row+1)*cellsize_),
            0            
          });
          g.push_back(arr3f{
            float(grid.minx_ + col*cellsize_),
            float(grid.miny_ + (row+1)*cellsize_),
            0            
          });
          auto lc = int(grid.getLinearCoord(row,col));
          tile_geom_cnts.push_back(tile_cnts[lc]);
          tile_geoms.push_back(g);
          tile_geom_ids.push_back(lc);
        }
      }
    };
  };

}