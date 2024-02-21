#pragma once

#include <geoflow/common.hpp>
#include <geoflow/geoflow.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <geoflow/parameters.hpp>

namespace geoflow::nodes::stepedge {
  typedef CGAL::Simple_cartesian<double>  K;
  typedef CGAL::Surface_mesh<K::Point_3>  SurfaceMesh;

  class MeshClipperNode:public Node {
    bool skip_clip_ = false;
    bool cgal_clip_ = false;
    bool smooth_normals_ = false;
    // float reduce_fraction_ = 0.5;
    // float agressiveness_ = 7.0;

    public:
    using Node::Node;
    void init() override {
      add_input("mesh", typeid(SurfaceMesh));
      add_input("bbox", typeid(Box));

      add_output("cgal_surface_mesh", typeid(SurfaceMesh));
      add_output("normals", typeid(vec3f));
      add_output("triangles", typeid(TriangleCollection));

      add_param(ParamBool(skip_clip_, "skip_clip", "Skip the clip"));
      add_param(ParamBool(cgal_clip_, "cgal_clip", "Use CGAL::Polygon_mesh_processing::clip instead of simpler but more robust triangle by triangle clip."));
      add_param(ParamBool(smooth_normals_, "smooth_normals", "Use use smooth vertex normals instead of flat face normals."));
      // add_param(ParamFloat(reduce_fraction_, "reduce_fraction", "Target reduction in nr of triangles"));
      // add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

  class MeshSimplifyNode:public Node {
    float stop_ratio_ = 0.5;
    bool border_correction_ = true;
    // float agressiveness_ = 7.0;

    public:
    using Node::Node;
    void init() override {
      add_input("cgal_surface_mesh", typeid(SurfaceMesh));

      add_output("cgal_surface_mesh", typeid(SurfaceMesh));

      add_param(ParamBool(border_correction_, "border_correction", "Correct ratio for border edges"));
      add_param(ParamBoundedFloat(stop_ratio_, 0, 1, "stop_ratio", "Target reduction ratio in nr of edges"));
      // add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

  class MeshSimplify2DNode:public Node {
    float error_ = 0.5;
    float minpts_ = 0.5;

    public:
    using Node::Node;
    void init() override {
      add_input("cgal_surface_mesh", typeid(SurfaceMesh));
      add_output("cgal_surface_mesh", typeid(SurfaceMesh));
      // add_output("wall_triangles", typeid(TriangleCollection));

      add_param(ParamBoundedFloat(error_, 0, 5, "error", "Target maximum eror after simplification"));
      add_param(ParamBoundedFloat(minpts_, 0, 10, "minpts", "Minimum number of elevation points per m2 inside a polygon"));
    }
    void process() override;
  };

  class MeshGridSimplifyNode:public Node {
    float cell_size_xy_ = 0.5;
    float cell_size_z_ = 0.5;
    // float agressiveness_ = 7.0;

    public:
    using Node::Node;
    void init() override {
      add_input("mesh", typeid(Mesh));
      add_input("bbox", typeid(Box));

      // add_output("tri", typeid(TriangleCollection));
      add_output("mesh", typeid(SurfaceMesh));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      add_param(ParamBoundedFloat(cell_size_xy_, 0, 1, "cell_size_xy", "Cellsize for x and y"));
      add_param(ParamBoundedFloat(cell_size_z_, 0, 1, "cell_size_z", "Cellsize for z"));
      // add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

  class Mesh2TriangleCollectionNode:public Node {
    // bool stop_ratio_;
    // float agressiveness_ = 7.0;

    public:
    using Node::Node;
    void init() override {
      add_input("cgal_surface_mesh", typeid(SurfaceMesh));

      add_output("triangles", typeid(TriangleCollection));
      add_output("normals", typeid(vec3f));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      // add_param(ParamBoundedFloat(stop_ratio_, 0, 1, "stop_ratio", "Target reduction ratio in nr of edges"));
      // add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };
  class Mesh2CGALSurfaceMeshNode:public Node {

    public:
    using Node::Node;
    void init() override {
      add_input("mesh", typeid(Mesh));
      add_output("cgal_surface_mesh", typeid(SurfaceMesh));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      // add_param(ParamBoundedFloat(stop_ratio_, 0, 1, "stop_ratio", "Target reduction ratio in nr of edges"));
      // add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };
  class SurfaceMesh2OFFNode:public Node {
    std::string filepath_="";

    public:
    using Node::Node;
    void init() override {
      add_input("cgal_surface_mesh", typeid(SurfaceMesh));

      add_param(ParamPath(filepath_, "filepath", "File path"));
    }
    void process() override;
  };

}