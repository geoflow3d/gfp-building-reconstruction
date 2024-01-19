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
#include "stepedge_nodes.hpp"
#include "MaxInscribedCircleNode.hpp"
#include "MeshSimplifyFastQuadNode.hpp"
#include "MeshProcessingNodes.hpp"

using namespace geoflow::nodes::stepedge;

void register_nodes(geoflow::NodeRegister& node_register) {
    node_register.register_node<AlphaShapeNode>("AlphaShape");
    node_register.register_node<PolygonExtruderNode>("PolygonExtruder");
    node_register.register_node<Arr2LinearRingsNode>("Arr2LinearRings");
    // node_register.register_node<ExtruderNode>("Extruder");
    node_register.register_node<LinearRingtoRingsNode>("LinearRingtoRings");
    // node_register.register_node<BuildArrFromRingsExactNode>("BuildArrFromRings");
    node_register.register_node<BuildArrFromLinesNode>("BuildArrFromLines");
    node_register.register_node<DetectLinesNode>("DetectLines");
    node_register.register_node<DetectPlanesNode>("DetectPlanes");
    node_register.register_node<LASInPolygonsNode>("LASInPolygons");
    #ifdef GFP_WITH_PDAL
      node_register.register_node<EptInPolygonsNode>("EptInPolygons");
    #endif
    node_register.register_node<BuildingSelectorNode>("BuildingSelector");
    node_register.register_node<ClusterLinesNode>("ClusterLines");
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
    node_register.register_node<ArrExtruderNode>("ArrExtruder");
    node_register.register_node<LOD1ExtruderNode>("LOD1Extruder");
    node_register.register_node<PolygonTriangulatorNode>("PolygonTriangulator");
    node_register.register_node<FacesSelectorNode>("FacesSelector");
    node_register.register_node<AttributeTesterNode>("AttributeTester");
    node_register.register_node<SkipReplaceTesterNode>("SkipReplaceTester");
    node_register.register_node<AttrRingsSelectorNode>("AttrRingsSelector");
    node_register.register_node<PolygonOffsetNode>("PolygonOffsetter");
    node_register.register_node<Arr2LinearRingsDebugNode>("Arr2LinearRingsDebug");
    node_register.register_node<DataCoverageCalcNode>("DataCoverageCalc");
    node_register.register_node<SnapRoundNode>("SnapRound");
    node_register.register_node<TriSnapNode>("TriSnap");
    node_register.register_node<BuildingRasteriseNode>("BuildingRasterise");
    node_register.register_node<PointCloudMergerNode>("PointCloudMerger");
    node_register.register_node<PC2PCDistancesCalculatorNode>("PC2PCDistancesCalculator");
    node_register.register_node<SegmentExtendNode>("SegmentExtend");
    node_register.register_node<PCFilterNode>("PCFilter");
    node_register.register_node<PlaneIntersectAllNode>("PlaneIntersectAll");
    node_register.register_node<RasterMergerNode>("RasterMerger");
    node_register.register_node<RoofPartitionRasteriseNode>("RoofPartitionRasterise");
    node_register.register_node<ClusterPointCloudNode>("ClusterPointCloud");
    node_register.register_node<ContourRegulariserNode>("ContourRegulariser");
    node_register.register_node<GridTilerNode>("GridTiler");
    node_register.register_node<MTC2MMNode>("MTC2MM");
    node_register.register_node<MaxInscribedCircleNode>("MaxInscribedCircle");
    node_register.register_node<MeshSimplifyFastQuadNode>("MeshSimplifyFastQuad");
    node_register.register_node<MeshSimplifyNode>("MeshSimplify");
    node_register.register_node<SurfaceMesh2OFFNode>("SurfaceMesh2OFF");
    node_register.register_node<MeshGridSimplifyNode>("MeshGridSimplify");
    node_register.register_node<MeshClipperNode>("MeshClipper");
    node_register.register_node<Mesh2TriangleCollectionNode>("Mesh2TriangleCollection");
    node_register.register_node<Mesh2CGALSurfaceMeshNode>("Mesh2CGALSurfaceMesh");
    node_register.register_node<RoofPartition3DBAGRasteriseNode>("RoofPartition3DBAGRasterise");
    node_register.register_node<MeshSimplify2DNode>("MeshSimplify2D");
}