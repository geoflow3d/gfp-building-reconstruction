cmake_minimum_required(VERSION 3.14)
project (building-reconstruction VERSION 0.0.1)

option(GFP_WITH_PDAL "Build parts dependent on PDAL" ON)

add_definitions(-DGF_PLUGIN_NAME=\"${CMAKE_PROJECT_NAME}\")

if(GFP_WITH_PDAL)
  add_definitions(-DGFP_WITH_PDAL)
endif(GFP_WITH_PDAL)

find_package(geoflow REQUIRED)

# ptinpoly
add_library(ptinpoly STATIC thirdparty/ptinpoly/ptinpoly.c)
set_target_properties(
  ptinpoly PROPERTIES 
  C_STANDARD 11
  POSITION_INDEPENDENT_CODE ON
)

if(EXISTS "${PROJECT_SOURCE_DIR}/.gitmodules")
  execute_process(
    COMMAND             git submodule update --init region-grower glm
    WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/thirdparty
    )
endif()

set(GLM_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/thirdparty/glm)

find_package(LASlib REQUIRED)
find_package(CGAL QUIET COMPONENTS Core REQUIRED)
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
  src/polygon_offset.cpp
  src/alpha_shape.cpp
  src/CityGMLMeshWriter.cpp
)

target_include_directories(gfp_buildingreconstruction PRIVATE
  src
  thirdparty/ptinpoly
  thirdparty/earcut
  thirdparty/region-grower
  ${PDAL_INCLUDE_DIRS}
  ${GLM_INCLUDE_DIRECTORIES}
)
target_link_libraries( gfp_buildingreconstruction PRIVATE
  geoflow-core
  ptinpoly
  CGAL::CGAL CGAL::CGAL_Core
  LASlib
  ${PDAL_LIBRARIES} 
)
target_link_directories( gfp_buildingreconstruction PRIVATE
  ${PDAL_LIBRARY_DIRS}
)
