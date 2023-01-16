#include <geoflow/geoflow.hpp>

namespace geoflow::nodes::stepedge {

  class MeshSimplifyFastQuadNode:public Node {
    float reduce_fraction_ = 0.5;
    float agressiveness_ = 7.0;

    public:
    using Node::Node;
    void init() override {
      add_input("triangles", typeid(TriangleCollection));

      add_output("triangles", typeid(TriangleCollection));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      add_param(ParamFloat(reduce_fraction_, "reduce_fraction", "Target reduction in nr of triangles"));
      add_param(ParamFloat(agressiveness_, "agressiveness", "Agressiveness"));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

}