cmake_minimum_required(VERSION 3.21)
project (building-reconstruction VERSION 0.4.5)

option(GFP_WITH_PDAL "Build parts dependent on PDAL" OFF)

add_definitions(-DGF_PLUGIN_NAME=\"${CMAKE_PROJECT_NAME}\")

if(PROJECT_IS_TOP_LEVEL)
  find_package(geoflow REQUIRED)
endif()

if(GFP_WITH_PDAL)
  add_definitions(-DGFP_WITH_PDAL)
endif(GFP_WITH_PDAL)

# ptinpoly
add_library(ptinpoly STATIC thirdparty/ptinpoly/ptinpoly.c)
set_target_properties(
  ptinpoly PROPERTIES 
  C_STANDARD 90
  POSITION_INDEPENDENT_CODE ON
)

if(EXISTS "${PROJECT_SOURCE_DIR}/.gitmodules")
  execute_process(
    COMMAND             git submodule update --init region-grower glm
    WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty
    )
endif()

set(GLM_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/thirdparty/glm)

find_package(laslib CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
add_definitions(-DCGAL_EIGEN3_ENABLED)
find_package(CGAL 5.4 QUIET COMPONENTS Core REQUIRED)
# PDAL
if(GFP_WITH_PDAL)
  find_package(PDAL REQUIRED)
endif(GFP_WITH_PDAL)

if (MSVC)
  # windows.h breaks std::min/std::max, fix by define
  add_definitions(-DNOMINMAX)
  # enable permissive compiling and/or statements
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive-")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:precise")
endif()

set(GF_PLUGIN_NAME ${PROJECT_NAME})
set(GF_PLUGIN_TARGET_NAME "gfp_buildingreconstruction")
set(GF_PLUGIN_REGISTER ${PROJECT_SOURCE_DIR}/register.hpp)
geoflow_create_plugin(
  src/stepedge_nodes.cpp
  src/node_optimise_arrangement.cpp
  src/points_in_polygons.cpp
  src/point_edge.cpp
  src/region_growing.cpp
  src/node_pcmesh_quality.cpp
  src/Raster.cpp
  src/heightfield_nodes.cpp
  src/polygon_triangulate.cpp
  src/line_regulariser.cpp
  src/build_arrangement_lines.cpp
  # src/build_arrangement_rings.cpp
  src/arrangement.cpp
  src/snap_round.cpp
  src/polygon_offset.cpp
  src/alpha_shape.cpp
  src/data_coverage_node.cpp
  src/detect_planes_node.cpp
  src/polygon_util.cpp
  src/cdt_util.cpp
  src/trisnap.cpp
  src/pip_util.cpp
  src/ClusterPointCloudNode.cpp
  src/ContourRegulariserNode.cpp
  src/MaxInscribedCircleNode.cpp
  src/MeshSimplifyFastQuadNode.cpp
  src/MeshClipperNode.cpp
  src/MeshSimplify.cpp
  src/MeshGridSimplify.cpp
  src/Mesh2TriangleCollectionNode.cpp
  src/tinsimp.cpp
)

target_include_directories(gfp_buildingreconstruction PRIVATE
  src
  thirdparty/ptinpoly
  thirdparty/earcut
  thirdparty/region-grower
  thirdparty/fast-quadric-mesh-simplification
  ${PDAL_INCLUDE_DIRS}
  ${GLM_INCLUDE_DIRECTORIES}
)
target_link_libraries( gfp_buildingreconstruction PRIVATE
  geoflow-core
  ptinpoly
  CGAL::CGAL CGAL::CGAL_Core Eigen3::Eigen
  LASlib
  ${PDAL_LIBRARIES} 
)
target_link_directories( gfp_buildingreconstruction PRIVATE
  ${PDAL_LIBRARY_DIRS}
)

if (MSVC)
  # collect dll's required for runtime
  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Release/ DESTINATION bin
    FILES_MATCHING 
    PATTERN "*.dll"
    PATTERN "gfp*" EXCLUDE)
endif()
