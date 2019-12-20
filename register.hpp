#include "stepedge_nodes.hpp"

using namespace geoflow::nodes::stepedge;

void register_nodes(geoflow::NodeRegister& node_register) {
    node_register.register_node<AlphaShapeNode>("AlphaShape");
    node_register.register_node<PolygonExtruderNode>("PolygonExtruder");
    node_register.register_node<Arr2LinearRingsNode>("Arr2LinearRings");
    node_register.register_node<ExtruderNode>("Extruder");
    node_register.register_node<LinearRingtoRingsNode>("LinearRingtoRings");
    node_register.register_node<BuildArrFromRingsExactNode>("BuildArrFromRings");
    node_register.register_node<BuildArrFromLinesNode>("BuildArrFromLines");
    node_register.register_node<DetectLinesNode>("DetectLines");
    node_register.register_node<DetectPlanesNode>("DetectPlanes");
    node_register.register_node<LASInPolygonsNode>("LASInPolygons");
    node_register.register_node<BuildingSelectorNode>("BuildingSelector");
    node_register.register_node<RegulariseLinesNode>("RegulariseLines");
    node_register.register_node<RegulariseRingsNode>("RegulariseRings");
    node_register.register_node<SimplifyPolygonNode>("SimplifyPolygon");
    node_register.register_node<Ring2SegmentsNode>("Ring2Segments");
    node_register.register_node<PrintResultNode>("PrintResult");
    node_register.register_node<PolygonGrowerNode>("PolygonGrower");
    node_register.register_node<PlaneIntersectNode>("PlaneIntersect");
    node_register.register_node<OptimiseArrangmentNode>("OptimiseArrangment");
    node_register.register_node<OptimiseArrangmentGridNode>("OptimiseArrangmentGrid");
    node_register.register_node<ArrDissolveNode>("ArrDissolve");
    node_register.register_node<PC2MeshQualityNode>("PC2MeshQuality");
    node_register.register_node<PCRasteriseNode>("PCRasterise");
    node_register.register_node<SegmentRasteriseNode>("SegmentRasterise");
    node_register.register_node<PolygonUnionNode>("PolygonUnion");
    node_register.register_node<Filter25DNode>("Filter25D");
    node_register.register_node<VecArr2LinearRingsNode>("VecArr2LinearRings");
}

namespace geoflow::nodes::stepedge {
  NodeRegisterHandle create_register() {
    auto R = NodeRegister::create(GF_PLUGIN_NAME);
    register_nodes(*R);
    return R;
  }
}