#include "stepedge_nodes.hpp"

namespace geoflow::nodes::stepedge {

  class LOD10GeneratorNode;
  class LOD13GeneratorNode;
  NodeRegisterHandle create_register() {
    auto R = NodeRegister::create("LoD13");
    R->register_node<AlphaShapeNode>("AlphaShape");
    R->register_node<PolygonExtruderNode>("PolygonExtruder");
    R->register_node<Arr2LinearRingsNode>("Arr2LinearRings");
    R->register_node<ExtruderNode>("Extruder");
    R->register_node<ProcessArrangementNode>("ProcessArrangement");
    R->register_node<LinearRingtoRingsNode>("LinearRingtoRings");
    R->register_node<BuildArrangementNode>("BuildArrangement");
    R->register_node<BuildArrFromRingsExactNode>("BuildArrFromRings");
    R->register_node<DetectLinesNode>("DetectLines");
    R->register_node<DetectPlanesNode>("DetectPlanes");
    R->register_node<ClassifyEdgePointsNode>("ClassifyEdgePoints");
    R->register_node<ComputeMetricsNode>("ComputeMetrics");
    R->register_node<LASInPolygonsNode>("LASInPolygons");
    R->register_node<BuildingSelectorNode>("BuildingSelector");
    R->register_node<RegulariseLinesNode>("RegulariseLines");
    R->register_node<RegulariseRingsNode>("RegulariseRings");
    R->register_node<SimplifyPolygonNode>("SimplifyPolygon");
    R->register_node<LOD10GeneratorNode>("LOD10Generator");
    R->register_node<LOD13GeneratorNode>("LOD13Generator");
    R->register_node<Ring2SegmentsNode>("Ring2Segments");
    R->register_node<PrintResultNode>("PrintResult");
    R->register_node<PolygonGrowerNode>("PolygonGrower");
    R->register_node<SegmentRasteriseNode>("SegmentRasterise");
    return R;
  }

  void create_lod13chart(NodeManager& N, bool use_linedetector) {
    auto R = create_register();

    auto DetectPlanes_node = N.create_node(R, "DetectPlanes", {300,75});
    N.name_node(DetectPlanes_node, "DetectPlanes_node");
    auto AlphaShape_node = N.create_node(R, "AlphaShape", {600,75});
    N.name_node(AlphaShape_node, "AlphaShape_node");
    auto DetectLines_node = N.create_node(R, "DetectLines", {900,-75});
    N.name_node(DetectLines_node, "DetectLines_node");
    auto RegulariseRings_node = N.create_node(R, "RegulariseRings", {1200,25});
    N.name_node(RegulariseRings_node, "RegulariseRings_node");
    auto BuildArrFromRings_node = N.create_node(R, "BuildArrFromRings", {1550,-125});
    N.name_node(BuildArrFromRings_node, "BuildArrFromRings_node");
    // auto ProcessArrangement_node = N.create_node(R, "ProcessArrangement");
    auto Arr2LinearRings_node = N.create_node(R, "Arr2LinearRings", {1550,225});
    N.name_node(Arr2LinearRings_node, "Arr2LinearRings_node");
    auto SimplifyPolygon_node = N.create_node(R, "SimplifyPolygon", {900,150});
    N.name_node(SimplifyPolygon_node, "SimplifyPolygon_node");
    // auto SimplifyPolygon_node_postfp = N.create_node(R, "SimplifyPolygon", {1200,-125});
    // N.name_node(SimplifyPolygon_node_postfp, "SimplifyPolygon_node_postfp");
    // auto SimplifyPolygon_node_postr = N.create_node(R, "SimplifyPolygon", {1200,175});
    // N.name_node(SimplifyPolygon_node_postr, "SimplifyPolygon_node_postr");
    auto Ring2Segments_node = N.create_node(R, "Ring2Segments", {900,50});
    N.name_node(Ring2Segments_node, "Ring2Segments_node");

    connect(DetectPlanes_node, AlphaShape_node, "pts_per_roofplane", "pts_per_roofplane");
    connect(DetectPlanes_node, BuildArrFromRings_node, "pts_per_roofplane", "pts_per_roofplane");
    
    if (!use_linedetector) {
      SimplifyPolygon_node->set_param("threshold_stop_cost", float(0.15));
      connect(AlphaShape_node, SimplifyPolygon_node, "alpha_rings", "polygons");
      connect(SimplifyPolygon_node, Ring2Segments_node, "polygons_simp", "rings");
      connect(Ring2Segments_node, RegulariseRings_node, "edge_segments", "edge_segments");
      connect(Ring2Segments_node, RegulariseRings_node, "ring_idx", "ring_idx");
      // connect(DetectLines_node, RegulariseRings_node, "edge_segments", "edge_segments");
      // connect(DetectLines_node, RegulariseRings_node, "ring_idx", "ring_idx");
      // connect(RegulariseRings_node, SimplifyPolygon_node_postr, "rings_out", "polygons");
      // connect(RegulariseRings_node, SimplifyPolygon_node_postfp, "footprint_out", "polygons");
      // connect(SimplifyPolygon_node_postr, BuildArrFromRings_node, "polygons_simp", "rings");
      // connect(SimplifyPolygon_node_postfp, BuildArrFromRings_node, "polygon_simp", "footprint");
    } else {
      connect(AlphaShape_node, DetectLines_node, "alpha_rings", "edge_points");
      connect(DetectLines_node, RegulariseRings_node, "edge_segments", "edge_segments");
      connect(DetectLines_node, RegulariseRings_node, "ring_idx", "ring_idx");
    }
    connect(RegulariseRings_node, BuildArrFromRings_node, "exact_rings_out", "rings");
    // if(regularise_footprint)
    connect(RegulariseRings_node, BuildArrFromRings_node, "exact_footprint_out", "footprint");
    connect(BuildArrFromRings_node, Arr2LinearRings_node, "arrangement", "arrangement");
  }

  class LOD13GeneratorNode:public Node {
    public:
    using Node::Node;
    void init() {
      add_input("point_clouds", typeid(std::vector<PointCollection>));
      add_input("polygons", typeid(LinearRingCollection));
      add_input_group("attributes", {typeid(vec1i), typeid(vec1f), typeid(vec1s), typeid(vec1b)});
      add_output("decomposed_footprints", typeid(LinearRingCollection));
//      add_output("attributes", typeid(AttributeMap));
      add_output("building_class", typeid(AttributeMap));
      add_output_group("attributes", {typeid(vec1i), typeid(vec1f), typeid(vec1s), typeid(vec1b)});

      add_param("step_height_threshold", (float) 2.0);
      // add_param("direct_alpharing", (bool) true);
      add_param("z_percentile", (float) 0.9);
      add_param("flood_to_unsegmented", (bool) true);
      add_param("dissolve_edges", (bool) true);
      add_param("dissolve_stepedges", (bool) true);
      add_param("use_only_hplanes", (bool) false);
      add_param("regularise_footprint", (bool) false);
      add_param("use_linedetector", (bool) false);
      // add_param("zrange_threshold", (float) 0.2);
      // add_param("merge_segid", (bool) true);
      // add_param("merge_zrange", (bool) false);
      // add_param("merge_step_height", (bool) true);
      // add_param("merge_unsegmented", (bool) false);
      // add_param("merge_dangling_egdes", (bool) false);
    }

    void gui(){
      // ImGui::Checkbox("direct_alpharing", &param<bool>("direct_alpharing"));
      
      ImGui::SliderFloat("Elevation percentile", &param<float>("z_percentile"), 0, 1);
      ImGui::Checkbox("Use only horizontal planes", &param<bool>("use_only_hplanes"));
      ImGui::Checkbox("Flood to unsegmented", &param<bool>("flood_to_unsegmented"));
      ImGui::Checkbox("Dissolve edges", &param<bool>("dissolve_edges"));
      ImGui::Checkbox("Dissolve stepedges", &param<bool>("dissolve_stepedges"));
      ImGui::SliderFloat("step_height_threshold", &param<float>("step_height_threshold"), 0, 100);
      ImGui::Checkbox("regularise_footprint", &param<bool>("regularise_footprint"));
      ImGui::Checkbox("use_linedetector", &param<bool>("use_linedetector"));
    }

    void process(){
      auto point_clouds = input("point_clouds").get<std::vector<PointCollection>>();
      auto polygons = input("polygons").get<LinearRingCollection>();
      
      
      // for each pair of polygon and point_cloud
        //create nodes and connections
        //run the thing
      if (point_clouds.size()!=polygons.size()) return;
      
      output_group("attributes").clear();
      
      for(auto& [name, term] : input_group("attributes").terminals) {
        auto& oT = output_group("attributes").add(name, term->connected_type);
        if(oT.type == typeid(vec1f))
          oT.set(vec1f());
        else if(oT.type == typeid(vec1i))
          oT.set(vec1i());
        else if(oT.type == typeid(vec1s))
          oT.set(vec1s());
        else if(oT.type == typeid(vec1b))
          oT.set(vec1b());
      }
      
      auto roof_pt_cnt_vec = vec1i();
      output_group("attributes").add("roof_pt_cnt_vec", typeid(vec1i));
      auto noid_a_vec = vec1f();
      output_group("attributes").add("noid_a", typeid(vec1f));
      auto noid_r_vec = vec1f();
      output_group("attributes").add("noid_r", typeid(vec1f));
      auto height_vec = vec1f();
      output_group("attributes").add("height", typeid(vec1f));
      auto bclass_vec = vec1i();
      output_group("attributes").add("bclass", typeid(vec1i));
      auto segid_vec = vec1i();
      output_group("attributes").add("segid", typeid(vec1i));

      LinearRingCollection all_cells;
      AttributeMap all_attributes;
      
      for(int i=0; i<point_clouds.size(); ++i) {
       std::cout << "13" << (param<bool>("use_only_hplanes") ? "H" : "A") << " fid: " << i << "\n";
        auto& points = point_clouds[i];
        auto& polygon = polygons[i];
        
        NodeManager N;
        create_lod13chart(N, param<bool>("use_linedetector"));

        // config and run
        // this should copy all parameters from this LOD13Generator node to the ProcessArrangement node
        N.nodes["BuildArrFromRings_node"]->set_params( dump_params(), true );
        N.nodes["BuildArrFromRings_node"]->set_param("extrude_unsegmented", param<bool>("use_only_hplanes"));
        N.nodes["RegulariseRings_node"]->set_param("regularise_fp", param<bool>("regularise_footprint"));
        N.nodes["DetectPlanes_node"]->set_param("only_horizontal", param<bool>("use_only_hplanes"));

        N.nodes["DetectPlanes_node"]->input("points").set(points);
        N.nodes["RegulariseRings_node"]->input("footprint").set(polygon);
        // if (!param<bool>("regularise_footprint")) {
          // N.nodes["BuildArrFromRings_node"]->input("footprint").set(polygon);
          // N.nodes["BuildArrFromRings_node"]->input("footprint").connected_type = TT_linear_ring;
        // }
        N.run(N.nodes["DetectPlanes_node"]);

        auto roof_pt_cnt = N.nodes["DetectPlanes_node"]->output("roof_pt_cnt").get<int>();
        auto classf = N.nodes["DetectPlanes_node"]->output("classf").get<float>();
        auto horiz = N.nodes["DetectPlanes_node"]->output("horiz_roofplane_cnt").get<float>();
        auto slant = N.nodes["DetectPlanes_node"]->output("slant_roofplane_cnt").get<float>();
        auto noseg_area_r = N.nodes["BuildArrFromRings_node"]->output("noseg_area_r").get<float>();
        auto noseg_area_a = N.nodes["BuildArrFromRings_node"]->output("noseg_area_a").get<float>();

        // note: the following will crash if the flowchart specified above is stopped halfway for some reason (eg missing output/connection)
        auto cells = N.nodes["Arr2LinearRings_node"]->output("linear_rings").get<LinearRingCollection>();
        auto attributes = N.nodes["Arr2LinearRings_node"]->output("attributes").get<AttributeMap>();

        for (int j=0; j<cells.size(); ++j) {
          all_cells.push_back(cells[j]);
          height_vec.push_back(attributes["height"][j]);
          segid_vec.push_back(attributes["segid"][j]);
          noid_a_vec.push_back(noseg_area_a);
          noid_r_vec.push_back(noseg_area_r);
          bclass_vec.push_back(classf);
          roof_pt_cnt_vec.push_back(roof_pt_cnt);
          for(auto& [name, iT] : input_group("attributes").terminals) {
            auto& oT = output_group("attributes").term(name);
            if(oT.type == typeid(vec1f)) {
              auto& veco = oT.get<vec1f&>();
              auto& veci = iT->get<vec1f&>();
              veco.push_back( veci[i] );
            }
            else if(oT.type == typeid(vec1i)) {
              auto& veco = oT.get<vec1i&>();
              auto& veci = iT->get<vec1i&>();
              veco.push_back( veci[i] );
            }
            else if(oT.type == typeid(vec1s)) {
              auto& veco = oT.get<vec1s&>();
              auto& veci = iT->get<vec1s&>();
              veco.push_back( veci[i] );
            }
          }
        }
      }
      output_group("attributes").term("roof_pt_cnt_vec").set(std::move(roof_pt_cnt_vec));
      output_group("attributes").term("noid_a").set(std::move(noid_a_vec));
      output_group("attributes").term("noid_r").set(std::move(noid_r_vec));
      output_group("attributes").term("height").set(std::move(height_vec));
      output_group("attributes").term("bclass").set(std::move(bclass_vec));
      output_group("attributes").term("segid").set(std::move(segid_vec));
      output("decomposed_footprints").set(std::move(all_cells));
    }
  };

  class LOD10GeneratorNode:public Node {
    public:
    using Node::Node;
    void init() {
      add_input("point_clouds", typeid(std::vector<PointCollection>));
      // add_output("attributes", typeid(AttributeMap));
      add_output_group("attributes", {typeid(vec1i), typeid(vec1f), typeid(vec1s), typeid(vec1b)});

      add_param("z_percentile", (float) 0.9);
    }

    void gui(){
      ImGui::SliderFloat("Elevation percentile", &param<float>("z_percentile"), 0, 1);
    }

    void process(){
      auto point_clouds = input("point_clouds").get<std::vector<PointCollection>>();

      // AttributeMap all_attributes;
      output_group("attributes").add("bclass", typeid(vec1i)).set(vec1i());
      auto& bclass_vec = output_group("attributes").term("bclass").get<vec1i&>();
      output_group("attributes").add("height", typeid(vec1f)).set(vec1f());
      auto& height_vec = output_group("attributes").term("height").get<vec1f&>();
      output_group("attributes").add("rms_error", typeid(vec1f)).set(vec1f());
      auto& rms_error_vec = output_group("attributes").term("rms_error").get<vec1f&>();
      output_group("attributes").add("roof_pt_cnt", typeid(vec1i)).set(vec1i());
      auto& roof_pt_cnt_vec = output_group("attributes").term("roof_pt_cnt").get<vec1i&>();

      auto R = create_register();
      
      for(int i=0; i<point_clouds.size(); i++) {
        std::cout << "10 fid: " << i << "\n";
        auto& point_cloud = point_clouds[i];
        
        NodeManager N;
        auto detect_planes_node = N.create_node(R, "DetectPlanes");
        detect_planes_node->input("points").set(point_cloud);

        N.run(detect_planes_node);

        auto classf = detect_planes_node->output("classf").get<float>();
        auto horiz = detect_planes_node->output("horiz_roofplane_cnt").get<float>();
        auto slant = detect_planes_node->output("slant_roofplane_cnt").get<float>();
        auto roof_pt_cnt = detect_planes_node->output("roof_pt_cnt").get<int>();
        bclass_vec.push_back(classf);
        roof_pt_cnt_vec.push_back(roof_pt_cnt);
        // all_attributes["horiz"].push_back(horiz);
        // all_attributes["slant"].push_back(slant);
        // note: the following will crash if the flowchart specified above is stopped halfway for some reason (eg missing output/connection)

        auto roofplanepts_per_fp = detect_planes_node->output("pts_per_roofplane").get<std::unordered_map<int, std::vector<Point>>>();

        std::vector<Point> points;
        for (auto& kv : roofplanepts_per_fp) {
          points.insert(points.end(), kv.second.begin(), kv.second.end());
        }
        if(points.size() < 2) {
          height_vec.push_back(0);
          rms_error_vec.push_back(0);
        } else {
          // if(polygons_feature.attr["height"][i]!=0) { //FIXME this is a hack!!
          std::sort(points.begin(), points.end(), [](linedect::Point& p1, linedect::Point& p2) {
            return p1.z() < p2.z();
          });
          auto elevation_id = int(param<float>("z_percentile")*float(points.size()-1));
          // std::cout << "id: " << elevation_id << ", size: " << points.size() << "\n";
          double elevation = points[elevation_id].z();
          double square_sum = 0;
          for (auto& p : points) {
            float d = elevation - p.z();
            square_sum += d*d;
          }
          rms_error_vec.push_back(CGAL::sqrt(square_sum/points.size()));
          height_vec.push_back(elevation);
        }

      }
      // output_group("attributes").term("bclass").set(std::move(bclass_vec));
      // output_group("attributes").term("height").set(std::move(height_vec));
      // output_group("attributes").term("rms_error").set(std::move(rms_error_vec));
    }
  };
}
