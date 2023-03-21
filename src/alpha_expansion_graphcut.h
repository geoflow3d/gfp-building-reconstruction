#ifndef CGAL_BOOST_GRAPH_ALPHA_EXPANSION_GRAPHCUT_H
// Copyright (c) 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Ilker O. Yaz, Simon Giraudot

#define CGAL_BOOST_GRAPH_ALPHA_EXPANSION_GRAPHCUT_H

#include <CGAL/Iterator_range.h>
#include <CGAL/assertions.h>
#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
#include <CGAL/Timer.h>
#endif

// #include <CGAL/boost/graph/Named_function_parameters.h>

#include <boost/version.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>

#if BOOST_VERSION >= 104400 // at this version kolmogorov_max_flow become depricated.
#  include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#else
#  include <boost/graph/kolmogorov_max_flow.hpp>
#endif

#include <vector>
#include <iostream>



namespace CGAL
{

/// \cond SKIP_IN_MANUAL
namespace internal
{

struct Alpha_expansion_old_API_wrapper_graph
{
  typedef std::size_t vertex_descriptor;
  typedef std::size_t edge_descriptor;
  typedef boost::directed_tag directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  typedef boost::edge_list_graph_tag traversal_category;
    
  typedef boost::counting_iterator<std::size_t> counting_iterator;
  typedef CGAL::Iterator_range<counting_iterator> counting_range;

  typedef CGAL::Identity_property_map<std::size_t> Vertex_index_map;
  typedef CGAL::Pointer_property_map<std::size_t>::type Vertex_label_map;

  struct Vertex_label_cost_map
  {
    typedef std::size_t key_type;
    typedef std::vector<double> value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;
      
    const std::vector<std::vector<double> >* cost_matrix;

    Vertex_label_cost_map (const std::vector<std::vector<double> >* cost_matrix)
      : cost_matrix (cost_matrix)
    { }

    friend reference get (const Vertex_label_cost_map& pmap, key_type idx)
    {
      std::vector<double> out;
      out.reserve (pmap.cost_matrix->size());
      for (std::size_t i = 0; i < pmap.cost_matrix->size(); ++ i)
        out.push_back ((*pmap.cost_matrix)[i][idx]);
      return out;
    }

  };

  typedef CGAL::Pointer_property_map<double>::const_type Edge_weight_map;
    
  const std::vector<std::pair<std::size_t, std::size_t> >& edges;
  const std::vector<double>& edge_weights;
  const std::vector<std::vector<double> >& cost_matrix;
  std::vector<std::size_t>& labels;

  Alpha_expansion_old_API_wrapper_graph (const std::vector<std::pair<std::size_t, std::size_t> >& edges,
                                         const std::vector<double>& edge_weights,
                                         const std::vector<std::vector<double> >& cost_matrix,
                                         std::vector<std::size_t>& labels)
    : edges (edges), edge_weights (edge_weights), cost_matrix (cost_matrix), labels (labels)
  { }

  friend counting_range vertices (const Alpha_expansion_old_API_wrapper_graph& graph)
  {
    return CGAL::make_range (boost::counting_iterator<std::size_t>(0),
                             boost::counting_iterator<std::size_t>(graph.labels.size()));
  }

  friend std::size_t num_vertices (const Alpha_expansion_old_API_wrapper_graph& graph) { return graph.labels.size(); }

  friend counting_range edges (const Alpha_expansion_old_API_wrapper_graph& graph)
  {
    return CGAL::make_range (boost::counting_iterator<std::size_t>(0),
                             boost::counting_iterator<std::size_t>(graph.edges.size()));
  }

  friend vertex_descriptor source (edge_descriptor ed, const Alpha_expansion_old_API_wrapper_graph& graph)
  { return graph.edges[ed].first; }
  friend vertex_descriptor target (edge_descriptor ed, const Alpha_expansion_old_API_wrapper_graph& graph)
  { return graph.edges[ed].second; }

  Vertex_index_map vertex_index_map() const { return Vertex_index_map(); }
  Vertex_label_map vertex_label_map() { return CGAL::make_property_map(labels); }
  Vertex_label_cost_map vertex_label_cost_map() const
  { return Vertex_label_cost_map(&cost_matrix); }
  Edge_weight_map edge_weight_map() const { return CGAL::make_property_map(edge_weights); }
};

////////////////////////////////////////////////////////////////////////////////////////
//   Comments about performance:
//
// 1) With BGL:
//     * Using adjacency_list:
//     ** Without pre-allocating vertex-list
//       | OutEdgeList | VertexList | Performance |
//       |    listS    |   listS    |     25.2    |
//       |    vecS     |   listS    |     22.7    |
//       |    listS    |   vecS     |     30.7    |
//       |    vecS     |   vecS     |     26.1    |
//
//     ** With pre-allocating vertex-list with max-node size
//        (Note: exact number of vertices are not certain at the beginning)
//       | OutEdgeList | VertexList | Performance |
//       |    listS    |   vecS     |     25.2    |
//       |    vecS     |   vecS     |     23.4    |
//
//     * Didn't try adjacency_matrix since our graph is sparse
//     ( Also one can check BGL book, performance section )
//
//    Decision:
//     * Alpha_expansion_graph_cut_boost: use adjacency_list<vecS, listS> without
//       pre-allocating vertex-list.
//
// 2) With Boykov-Kolmogorov MAXFLOW software:
//   (http://pub.ist.ac.at/~vnk/software/maxflow-v2.21.src.tar.gz)
//                                  | Performance |
//                                  |     3.1     |
//     * Alpha_expansion_graph_cut_boykov_kolmogorov provides an implementation.
//       MAXFLOW does not provide any option for pre-allocation (It is possible with v_3.02 though).
//
// Typical Benchmark result provided by Ilker
//                                 | construction of vertices  |  construction of edges    | graph cut  | Total
//   -----------------------------------------------------------------------------------------------------------
//   boost with an adjacency list  |         1.53              | 1.51                      |  3.00      | 6.04
//   boost with CSR                | 0.11 (gather in a vector) | 0.15 (gather in a vector) |  2.67      | 2.93
//   MaxFlow                       |       0.042               | 0.076                     |  1.043     | 1.161
//
// The main issue for now with CSR is the construction of the opposite edge map that is too costly,
// since it is done by exploring all edges to find opposite
////////////////////////////////////////////////////////////////////////////////////////

} // namespace internal

/**
 * @brief Implements alpha-expansion graph cut algorithm.
 *
 * For representing graph, it uses adjacency_list with OutEdgeList = vecS, VertexList = listS.
 * Also no pre-allocation is made for vertex-list.
 */
class Alpha_expansion_boost_adjacency_list_tag
{
private:
  typedef boost::adjacency_list_traits<boost::vecS, boost::listS, boost::directedS>
  Adjacency_list_traits;

  typedef boost::adjacency_list<boost::vecS, boost::listS, boost::directedS,
          // 4 vertex properties
          boost::property<boost::vertex_index_t, std::size_t,
          boost::property<boost::vertex_color_t, boost::default_color_type,
          boost::property<boost::vertex_distance_t, double,
          boost::property<boost::vertex_predecessor_t, Adjacency_list_traits::edge_descriptor >
          > > >,
          // 3 edge properties
          boost::property<boost::edge_capacity_t, double,
          boost::property<boost::edge_residual_capacity_t, double,
          boost::property<boost::edge_reverse_t, Adjacency_list_traits::edge_descriptor> >
          > > Graph;

  typedef boost::graph_traits<Graph> Traits;
  typedef boost::color_traits<boost::default_color_type> ColorTraits;

public:

  typedef Traits::vertex_descriptor Vertex_descriptor;
  typedef Traits::vertex_iterator   Vertex_iterator;
  typedef Traits::edge_descriptor   Edge_descriptor;
  typedef Traits::edge_iterator     Edge_iterator;

private:

  Graph graph;
  Vertex_descriptor cluster_source;
  Vertex_descriptor cluster_sink;
  
public:

  void clear_graph()
  {
    graph.clear();
    cluster_source = boost::add_vertex(graph);
    cluster_sink = boost::add_vertex(graph);
  }

  Vertex_descriptor add_vertex()
  {
    return boost::add_vertex(graph);
  }

  void add_tweight (Vertex_descriptor& v, double w1, double w2)
  {
    add_edge (cluster_source, v, w1, 0);
    add_edge (v, cluster_sink, w2, 0);
  }

  void init_vertices()
  {
    // initialize vertex indices, it is necessary since we are using VertexList = listS
    Vertex_iterator v_begin, v_end;
    Traits::vertices_size_type index = 0;
    for(boost::tie(v_begin, v_end) = vertices(graph); v_begin != v_end; ++v_begin) {
      boost::put(boost::vertex_index, graph, *v_begin, index++);
    }
  }

  double max_flow()
  {
#if BOOST_VERSION >= 104400
    return boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                                                      cluster_sink);
#else
    return boost::kolmogorov_max_flow(graph, cluster_source, cluster_sink);
#endif
  }

  template <typename VertexLabelMap, typename InputVertexDescriptor>
  void update(VertexLabelMap vertex_label_map,
              const std::vector<Vertex_descriptor>& inserted_vertices,
              InputVertexDescriptor vd,
              std::size_t vertex_i,
              std::size_t alpha)
  {
    boost::default_color_type color = boost::get(boost::vertex_color, graph,
                                                 inserted_vertices[vertex_i]);
    if(get (vertex_label_map, vd) != alpha
       && color == ColorTraits::white()) //new comers (expansion occurs)
      put (vertex_label_map, vd, alpha);
  }
  
  void add_edge (Vertex_descriptor& v1, Vertex_descriptor& v2, double w1, double w2)
  {
    Edge_descriptor v1_v2, v2_v1;
    bool v1_v2_added, v2_v1_added;

    boost::tie(v1_v2, v1_v2_added) = boost::add_edge(v1, v2, graph);
    boost::tie(v2_v1, v2_v1_added) = boost::add_edge(v2, v1, graph);

    CGAL_assertion(v1_v2_added && v2_v1_added);
    //put edge capacities
    boost::put(boost::edge_reverse, graph, v1_v2, v2_v1);
    boost::put(boost::edge_reverse, graph, v2_v1, v1_v2);

    //map reverse edges
    boost::put(boost::edge_capacity, graph, v1_v2, w1);
    boost::put(boost::edge_capacity, graph, v2_v1, w2);

  }
};

// another implementation using compressed_sparse_row_graph
// for now there is a performance problem while setting reverse edges
// if that can be solved, it is faster than Alpha_expansion_graph_cut_boost
class Alpha_expansion_boost_compressed_sparse_row_tag
{
private:
  // CSR only accepts bundled props
  struct VertexP {
    boost::default_color_type vertex_color;
    double vertex_distance_t;
    // ? do not now there is another way to take it, I think since edge_descriptor does not rely on properties
    // this should be fine...
    boost::compressed_sparse_row_graph<boost::directedS>::edge_descriptor
    vertex_predecessor;
  };

  struct EdgeP {
    double edge_capacity;
    double edge_residual_capacity;
    boost::compressed_sparse_row_graph<boost::directedS>::edge_descriptor
    edge_reverse;
  };

  typedef boost::compressed_sparse_row_graph<boost::directedS,
          VertexP, EdgeP> Graph;

  typedef boost::graph_traits<Graph> Traits;
  typedef boost::color_traits<boost::default_color_type> ColorTraits;

public:
  
  typedef Traits::vertex_descriptor Vertex_descriptor;
  typedef Traits::vertex_iterator   Vertex_iterator;
  typedef Traits::edge_descriptor   Edge_descriptor;
  typedef Traits::edge_iterator     Edge_iterator;

private:

  Graph graph;
  std::size_t nb_vertices;
  std::vector<std::pair<std::size_t, std::size_t> > edge_map;
  std::vector<EdgeP>                edge_map_weights;

public:
  void clear_graph()
  {
    nb_vertices = 2;
    edge_map.clear();
    edge_map_weights.clear();
    // edge_map.reserve(labels.size() *
    //                  8); // there is no way to know exact edge count, it is a heuristic value
    // edge_map_weights.reserve(labels.size() * 8);
  }

  Vertex_descriptor add_vertex()
  {
    return (nb_vertices ++);
  }

  void add_tweight (Vertex_descriptor& v, double w1, double w2)
  {
    add_edge (0, v, w1, 0);
    add_edge (v, 1, w2, 0);
  }

  void init_vertices()
  {
#if BOOST_VERSION >= 104000
    graph = Graph(boost::edges_are_unsorted, edge_map.begin(), edge_map.end(),
                  edge_map_weights.begin(), nb_vertices);
#else
    graph= Graph(edge_map.begin(), edge_map.end(),
                 edge_map_weights.begin(), nb_vertices);
#endif

    // PERFORMANCE PROBLEM
    // need to set reverse edge map, I guess there is no way to do that before creating the graph
    // since we do not have edge_descs
    // however from our edge_map, we know that each (2i, 2i + 1) is reverse pairs, how to facilitate that ?
    // will look it back
    Graph::edge_iterator ei, ee;
    for(boost::tie(ei, ee) = boost::edges(graph); ei != ee; ++ei) {
      Graph::vertex_descriptor v1 = boost::source(*ei, graph);
      Graph::vertex_descriptor v2 = boost::target(*ei, graph);
      std::pair<Graph::edge_descriptor, bool> opp_edge = boost::edge(v2, v1, graph);

      CGAL_assertion(opp_edge.second);
      graph[opp_edge.first].edge_reverse =
        *ei; // and edge_reverse of *ei will be (or already have been) set by the opp_edge
    }
  }

  double max_flow()
  {
#if BOOST_VERSION >= 104400
    // since properties are bundled, defaults does not work need to specify them
    return boost::boykov_kolmogorov_max_flow
      (graph,
       boost::get(&EdgeP::edge_capacity, graph),
       boost::get(&EdgeP::edge_residual_capacity, graph),
       boost::get(&EdgeP::edge_reverse, graph),
       boost::get(&VertexP::vertex_predecessor, graph),
       boost::get(&VertexP::vertex_color, graph),
       boost::get(&VertexP::vertex_distance_t, graph),
       boost::get(boost::vertex_index,
                  graph), // this is not bundled, get it from graph (CRS provides one)
       0, 1);
#else
    return boost::kolmogorov_max_flow
       (graph,
        boost::get(&EdgeP::edge_capacity, graph),
        boost::get(&EdgeP::edge_residual_capacity, graph),
        boost::get(&EdgeP::edge_reverse, graph),
        boost::get(&VertexP::vertex_predecessor, graph),
        boost::get(&VertexP::vertex_color, graph),
        boost::get(&VertexP::vertex_distance_t, graph),
        boost::get(boost::vertex_index,
                   graph), // this is not bundled, get it from graph
        0, 1);
#endif
  }

  template <typename VertexLabelMap, typename InputVertexDescriptor>
  void update(VertexLabelMap vertex_label_map,
              const std::vector<Vertex_descriptor>& inserted_vertices,
              InputVertexDescriptor vd,
              std::size_t vertex_i,
              std::size_t alpha)
  {
    boost::default_color_type color =  graph[vertex_i + 2].vertex_color;
    if(get(vertex_label_map, vd)!= alpha
       && color == ColorTraits::white()) //new comers (expansion occurs)
      put(vertex_label_map, vd, alpha);
  }
  
  void add_edge(Vertex_descriptor v1, Vertex_descriptor v2, double w1, double w2)
  {
    edge_map.push_back(std::make_pair(v1, v2));
    EdgeP p1;
    p1.edge_capacity = w1;
    edge_map_weights.push_back(p1);

    edge_map.push_back(std::make_pair(v2, v1));
    EdgeP p2;
    p2.edge_capacity = w2;
    edge_map_weights.push_back(p2);
  }
  

};

// Default version using boost adjacency list
template <typename InputGraph,
          typename EdgeWeightMap,
          typename VertexIndexMap,
          typename VertexLabelCostMap,
          typename VertexLabelMap>
double alpha_expansion_graphcut (const InputGraph& input_graph,
                                 EdgeWeightMap edge_weight_map,
                                 VertexIndexMap vertex_index_map,
                                 VertexLabelCostMap vertex_label_cost_map,
                                 VertexLabelMap vertex_label_map)
{
  return alpha_expansion_graphcut
    (input_graph, edge_weight_map, vertex_index_map, vertex_label_cost_map, vertex_label_map,
     Alpha_expansion_boost_adjacency_list_tag());
}

/// \endcond

// NOTE: latest performances check (2019-07-22)
//
// Using a random graph with 50000 vertices, 100000 edges and 30 labels:
//
// METHOD                 TIMING           MEMORY
// Boost Adjacency list   49s              122MiB
// Boost CSR              187s             77MiB
// MaxFlow                12s              717MiB

/**
   \ingroup PkgBGLPartition

   Regularizes a partition of a graph into `n` labels using the Alpha
   Expansion algorithm \cgalCite{Boykov2001FastApproximate}.

   For a graph \f$(V,E)\f$, this function seeks for the partioning `f`
   that minimizes the following energy:

   \f[
   \mathrm{E}(f) = \sum_{\{v0,v1\} \in E} W(v0,v1) + \sum_{v \in V} C(f_v)
   \f]

   where \f$W(v0,v1)\f$ is the weight associated to the edge
   \f$\{v0,v1\}\f$ and \f$C(f_v)\f$ is the cost of assigning the
   vertex \f$v\f$ to the labeling \f$f\f$.

   \tparam InputGraph a model of `Graph`

   \tparam EdgeWeightMap a model of `ReadablePropertyMap` with
   `boost::graph_traits<InputGraph>::%edge_descriptor` as key and `double`
   as value

   \tparam VertexIndexMap a model of `ReadablePropertyMap` with
   `boost::graph_traits<InputGraph>::%vertex_descriptor` as key and
   `std::size_t` as value

   \tparam VertexLabelCostMap a model of `ReadablePropertyMap`
   with `boost::graph_traits<InputGraph>::%vertex_descriptor` as key and
   `std::vector<double>` as value

   \tparam VertexLabelMap a model of `ReadWritePropertyMap` with 
   `boost::graph_traits<InputGraph>::%vertex_descriptor` as key and
   `std::size_t` as value

   \tparam AlphaExpansionImplementationTag optional tag used to select
   which implementation of the Alpha Expansion should be
   used. Available implementation tags are:

   - `CGAL::Alpha_expansion_boost_adjacency_list` (default)
   - `CGAL::Alpha_expansion_boost_compressed_sparse_row_tag`
   - `CGAL::Alpha_expansion_MaxFlow_tag`

   \note The `MaxFlow` implementation is provided by the \ref PkgSurfaceMeshSegmentationRef
   under GPL license. The header
   `<CGAL/boost/graph/Alpha_expansion_MaxFlow_tag.h>`
   must be included if users want to use this implementation.

   \param input_graph the input graph.

   \param edge_weight_map a property map providing the weight of each
   edge.

   \param vertex_index_map a property map providing the index of each
   vertex.

   \param vertex_label_map a property map providing the label of each
   vertex. This map will be updated by the algorithm with the
   regularized version of the partitioning.

   \param vertex_label_cost_map a property map providing, for each
   vertex, an `std::vector` containing the cost of this vertex to
   belong to each label. For example,
   `get(vertex_label_cost_map, vd)[label_idx]` returns the cost
   of vertex `vd` to belong to the label `label_idx`.
*/
template <typename InputGraph,
          typename EdgeWeightMap,
          typename VertexIndexMap,
          typename VertexLabelCostMap,
          typename VertexLabelMap,
          typename AlphaExpansionImplementationTag>
double alpha_expansion_graphcut (const InputGraph& input_graph,
                                 EdgeWeightMap edge_weight_map,
                                 VertexIndexMap vertex_index_map,
                                 VertexLabelCostMap vertex_label_cost_map,
                                 VertexLabelMap vertex_label_map,
                                 const AlphaExpansionImplementationTag&,
                                 size_t n_iterations=5)
{
  typedef boost::graph_traits<InputGraph> GT;
  typedef typename GT::edge_descriptor input_edge_descriptor;
  typedef typename GT::vertex_descriptor input_vertex_descriptor;

  typedef AlphaExpansionImplementationTag Alpha_expansion;
  typedef typename Alpha_expansion::Vertex_descriptor Vertex_descriptor;

  Alpha_expansion alpha_expansion;

  // TODO: check this hardcoded parameter
  const double tolerance = 1e-10;

  double min_cut = (std::numeric_limits<double>::max)();

#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
  double vertex_creation_time, edge_creation_time, cut_time;
  vertex_creation_time = edge_creation_time = cut_time = 0.0;
#endif

  std::vector<Vertex_descriptor> inserted_vertices;
  inserted_vertices.resize(num_vertices (input_graph));

  std::size_t number_of_labels = get(vertex_label_cost_map, *vertices(input_graph).first).size();

  bool success;
  size_t iteration = 0;
  do {
    success = false;

    for (std::size_t alpha = 0; alpha < number_of_labels; ++ alpha)
    {
      alpha_expansion.clear_graph();

#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
      Timer timer;
      timer.start();
#endif

      // For E-Data
      // add every input vertex as a vertex to the graph, put edges to source & sink vertices
      auto it_pair = vertices(input_graph);
      for (auto vd_ = it_pair.first; vd_ != it_pair.second; ++vd_ )
      {
        input_vertex_descriptor vd = *vd_;
        std::size_t vertex_i = get(vertex_index_map, vd);
        Vertex_descriptor new_vertex = alpha_expansion.add_vertex();
        inserted_vertices[vertex_i] = new_vertex;
        double source_weight = get(vertex_label_cost_map, vd)[alpha];
        // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
        // making sink_weight 'infinity' guarantee this.
        double sink_weight = (get(vertex_label_map, vd) == alpha ?
                              (std::numeric_limits<double>::max)()
                              : get(vertex_label_cost_map, vd)[get(vertex_label_map, vd)]);

        alpha_expansion.add_tweight(new_vertex, source_weight, sink_weight);
      }
#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
      vertex_creation_time += timer.time();
      timer.reset();
#endif

      // For E-Smooth
      // add edge between every vertex,
      auto it_pair_edges = edges(input_graph);
      for (auto ed_ = it_pair_edges.first; ed_ != it_pair_edges.second; ++ed_)
      {
        input_edge_descriptor ed = *ed_;
        input_vertex_descriptor vd1 = source(ed, input_graph);
        input_vertex_descriptor vd2 = target(ed, input_graph);
        std::size_t idx1 = get (vertex_index_map, vd1);
        std::size_t idx2 = get (vertex_index_map, vd2);
          
        double weight = get (edge_weight_map, ed);
          
        Vertex_descriptor v1 = inserted_vertices[idx1],
          v2 = inserted_vertices[idx2];
          
        std::size_t label_1 = get (vertex_label_map, vd1);
        std::size_t label_2 = get (vertex_label_map, vd2);
        if(label_1 == label_2) {
          if(label_1 != alpha) {
            alpha_expansion.add_edge(v1, v2, weight, weight);
          }
        } else {
          Vertex_descriptor inbetween = alpha_expansion.add_vertex();

          double w1 = (label_1 == alpha) ? 0 : weight;
          double w2 = (label_2 == alpha) ? 0 : weight;
          alpha_expansion.add_edge(inbetween, v1, w1, w1);
          alpha_expansion.add_edge(inbetween, v2, w2, w2);
          alpha_expansion.add_tweight(inbetween, 0., weight);
        }
      }
#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
      edge_creation_time += timer.time();
#endif

      alpha_expansion.init_vertices();

#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
      timer.reset();
#endif
      
      double flow = alpha_expansion.max_flow();
      
#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
      cut_time += timer.time();
#endif

      // std::cerr<< "min_cut: " << min_cut << ", flow: " << flow << " [a=" << alpha << "]\n";
      if(min_cut - flow <= flow * tolerance) {
        continue;
      }
      min_cut = flow;
      success = true;
      //update labeling
      auto it_pair_v = vertices(input_graph);
      for ( auto vd_ = it_pair_v.first; vd_ != it_pair_v.second; ++vd_)
      {
        input_vertex_descriptor vd = *vd_;
        std::size_t vertex_i = get (vertex_index_map, vd);
        alpha_expansion.update(vertex_label_map, inserted_vertices, vd, vertex_i, alpha);
      }
    }
  // } while(success);
  } while(++iteration < n_iterations);

#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
  CGAL_TRACE_STREAM << "vertex creation time: " << vertex_creation_time <<
    std::endl;
  CGAL_TRACE_STREAM << "edge creation time: " << edge_creation_time << std::endl;
  CGAL_TRACE_STREAM << "max flow algorithm time: " << cut_time << std::endl;
#endif

  return min_cut;
}


/// \cond SKIP_IN_MANUAL
// Old API
inline double alpha_expansion_graphcut (const std::vector<std::pair<std::size_t, std::size_t> >& edges,
                                        const std::vector<double>& edge_weights,
                                        const std::vector<std::vector<double> >& cost_matrix,
                                        std::vector<std::size_t>& labels)
{
  internal::Alpha_expansion_old_API_wrapper_graph graph (edges, edge_weights, cost_matrix, labels);
  
  return alpha_expansion_graphcut(graph,
                                  graph.edge_weight_map(),
                                  graph.vertex_index_map(),
                                  graph.vertex_label_cost_map(),
                                  graph.vertex_label_map());
}

template <typename AlphaExpansionImplementationTag>
double alpha_expansion_graphcut (const std::vector<std::pair<std::size_t, std::size_t> >& edges,
                                 const std::vector<double>& edge_weights,
                                 const std::vector<std::vector<double> >& cost_matrix,
                                 std::vector<std::size_t>& labels,
                                 const AlphaExpansionImplementationTag&)
{
  internal::Alpha_expansion_old_API_wrapper_graph graph (edges, edge_weights, cost_matrix, labels);
  
  return alpha_expansion_graphcut(graph,
                                  graph.edge_weight_map(),
                                  graph.vertex_index_map(),
                                  graph.vertex_label_cost_map(),
                                  graph.vertex_label_map(),
                                  AlphaExpansionImplementationTag());
}
/// \endcond


}//namespace CGAL
#endif //CGAL_BOOST_GRAPH_ALPHA_EXPANSION_GRAPHCUT_H
