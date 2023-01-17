#pragma once

#include <geoflow/geoflow.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

namespace geoflow::nodes::stepedge {
  typedef CGAL::Simple_cartesian<double>  K;
  typedef CGAL::Surface_mesh<K::Point_3>  SurfaceMesh;

  class MeshClipperNode:public Node {
    // float reduce_fraction_ = 0.5;
    // float agressiveness_ = 7.0;

    public:
    using Node::Node;
    void init() override {
      add_input("mesh", typeid(Mesh));
      add_input("bbox", typeid(Box));

      add_output("triangles", typeid(TriangleCollection));
      add_output("cgal_surface_mesh", typeid(SurfaceMesh));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      // add_param(ParamFloat(reduce_fraction_, "reduce_fraction", "Target reduction in nr of triangles"));
      // add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

  class MeshSimplifyNode:public Node {
    float stop_ratio_ = 0.5;
    // float agressiveness_ = 7.0;

    public:
    using Node::Node;
    void init() override {
      add_input("cgal_surface_mesh", typeid(SurfaceMesh));

      add_output("triangles", typeid(TriangleCollection));
      add_output("cgal_surface_mesh", typeid(SurfaceMesh));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      add_param(ParamBoundedFloat(stop_ratio_, 0, 1, "stop_ratio", "Target reduction ratio in nr of edges"));
      // add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

}