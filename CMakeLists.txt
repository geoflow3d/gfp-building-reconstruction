cmake_minimum_required(VERSION 3.14)
project (building-reconstruction VERSION 0.0.1)

add_definitions(-DGF_PLUGIN_NAME=\"${CMAKE_PROJECT_NAME}\")

find_package(geoflow REQUIRED)

# ptinpoly
add_library(ptinpoly STATIC thirdparty/ptinpoly/ptinpoly.c)
set_target_properties(
  ptinpoly PROPERTIES 
  C_STANDARD 11
  POSITION_INDEPENDENT_CODE ON
)

find_package(glm)
find_package(LASlib)
find_package(CGAL QUIET COMPONENTS Core)
#earcut
include_directories(
  src 
  thirdparty/ptinpoly
  thirdparty/earcut
  thirdparty/region-grower
)

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:precise")
endif()

add_library(gfp_stepedge MODULE
  plugin.cpp
  src/stepedge_nodes.cpp
  src/node_optimise_arrangement.cpp
  src/point_edge.cpp
  src/region_growing.cpp
  src/node_pcmesh_quality.cpp
  src/Raster.cpp
  src/heightfield_nodes.cpp
)
target_link_libraries( gfp_stepedge
  geoflow-core
  ptinpoly
  glm
  CGAL::CGAL CGAL::CGAL_Core
  LASlib
)
set_target_properties(
  gfp_stepedge
  PROPERTIES CXX_STANDARD 17
)
install(TARGETS gfp_stepedge DESTINATION ${geoflow_PLUGIN_FOLDER})
