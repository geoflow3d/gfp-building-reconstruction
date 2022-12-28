#include <geoflow/geoflow.hpp>

namespace geoflow::nodes::stepedge {

  class MaxInscribedCircleNode:public Node {
    float cellsize_ = 2.0;
    float polygon_densify_ = 0.5;

    public:
    using Node::Node;
    void init() override {
      add_input("polygon", typeid(LinearRing));
      add_input("pointcloud", typeid(PointCollection));

      add_output("max_diameter", typeid(float));
      add_output("max_circle", typeid(LinearRing));
      add_output("vd_pts", typeid(PointCollection));

      // add_param(ParamBool(flatten, "flatten", "Ignore Z coordinates in clustering"));
      add_param(ParamFloat(cellsize_, "cellsize", "Cellsize"));
      add_param(ParamFloat(polygon_densify_, "polygon_densify", "Densify polygon edges using this threshold."));
      // add_param(ParamInt(metrics_normal_k, "metrics_normal_k", "Number of neighbours used for normal estimation"));
    }
    void process() override;
  };

}