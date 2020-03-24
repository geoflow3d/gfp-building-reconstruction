#include "stepedge_nodes.hpp"

#include <sstream>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/Options.hpp>
#include <pdal/io/EptReader.hpp>
#include <pdal/filters/RangeFilter.hpp>
#include <pdal/filters/CropFilter.hpp>

#include <lasreader.hpp>

#include "ptinpoly.h"

namespace geoflow::nodes::stepedge {

pGridSet build_grid(const vec3f& ring) {
  int Grid_Resolution = 20;

  int size = ring.size();
  std::vector<pPipoint> pgon;
  for (auto& p : ring) {
    pgon.push_back(new Pipoint{ p[0],p[1] });
  }
  pGridSet grid_set = new GridSet();
  // skip last point in the ring, ie the repetition of the first vertex
  GridSetup(&pgon[0], pgon.size(), Grid_Resolution, grid_set);
  for (int i = 0; i < size; i++) {
    delete pgon[i];
  }
  return grid_set;
}

class PointsInPolygonsCollector  {
  gfSingleFeatureInputTerminal& polygons;
  gfSingleFeatureOutputTerminal& point_clouds;
  gfSingleFeatureOutputTerminal& ground_elevations;

  // ground elevations
  std::vector<std::vector<float>> point_clouds_ground;
  RasterTools::Raster pindex;
  std::vector<std::vector<size_t>> pindex_vals;
  std::vector<pGridSet> poly_grids;
  
  public:
  Box completearea_bb;
  float min_ground_elevation = std::numeric_limits<float>::max();

  PointsInPolygonsCollector(gfSingleFeatureInputTerminal& polygons, gfSingleFeatureOutputTerminal& point_clouds, gfSingleFeatureOutputTerminal& ground_elevations, float& cellsize, float& buffer)
  : polygons(polygons), point_clouds(point_clouds), ground_elevations(ground_elevations) {
    point_clouds_ground.resize(polygons.size());
    point_clouds.resize<PointCollection>(polygons.size());

    // make a vector of BOX2D for the set of input polygons
    // build point in polygon grids
    for (size_t i=0; i<polygons.size(); ++i) {
      auto ring = polygons.get<LinearRing>(i);
      poly_grids.push_back(build_grid(ring));
      completearea_bb.add(ring.box());
    }

    // build an index grid for the polygons

    // build raster index (pindex) and store for each raster cell all the polygons with intersecting bbox (pindex_vals)
    float minx = completearea_bb.min()[0]-buffer; 
    float miny = completearea_bb.min()[1]-buffer;
    float maxx = completearea_bb.max()[0]+buffer;
    float maxy = completearea_bb.max()[1]+buffer;
    pindex =  RasterTools::Raster(cellsize, minx, maxx, miny, maxy);
    pindex_vals.resize(pindex.dimx_*pindex.dimy_);

    // populate pindex_vals
    for (size_t i=0; i<polygons.size(); ++i) {
      auto ring = polygons.get<LinearRing>(i);
      auto& b = ring.box();
      size_t r_min = pindex.getRow(b.min()[0], b.min()[1]);
      size_t c_min = pindex.getCol(b.min()[0], b.min()[1]);
      size_t r_max = static_cast<size_t>( ceil(b.size_y() / cellsize) ) + r_min;
      size_t c_max = static_cast<size_t>( ceil(b.size_x() / cellsize) ) + c_min;
      if (r_max >= pindex.dimy_)
        r_max = pindex.dimy_-1;
      if (c_max >= pindex.dimx_)
        c_max = pindex.dimx_-1;
      for(size_t r = r_min; r <= r_max; ++r) {
        for(size_t c = c_min; c <= c_max; ++c) {
          pindex_vals[r*pindex.dimx_ + c].push_back(i);
        }
      }
    }

  }

  ~PointsInPolygonsCollector() {
    for (int i=0; i<poly_grids.size(); i++) {
      delete poly_grids[i];
    }
  }

  void add_point(arr3f point, int point_class) {
    // look up grid index cell and do pip for all polygons retreived from that cell
    size_t lincoord = pindex.getLinearCoord(point[0],point[1]);
    if (lincoord >= pindex_vals.size() || lincoord < 0) {
      // std::cout << "Point (" << point[0] << ", " <<point[1] << ", "  << point[2] << ") is not in the polygon bbox.\n";
      return;
    }
    // For the ground points we only test if the point is within the grid
    // index cell, but we do not do a pip for the footprint itself.
    // The reasons for this:
    //  - If we would do a pip for with the ground points, we would need to
    //    have a separate list of buffered footprints that we use for the
    //    ground pip. Still it would be often the case that there are no
    //    ground points found for a buffered footprint. Thus we need a larger
    //    area in which we can guarantee that we find at least a couple of
    //    ground points.
    //  - Because we work with Netherlands data, the ground relief is small.
    //    Thus a single ground height value per grid cell is good enough for
    //    representing the ground/floor elevation of the buildings in that
    //    grid cell.
    pPipoint pipoint = new Pipoint{point[0],point[1]};
    for(size_t& poly_i : pindex_vals[lincoord]) {
      if (point_class == 2) {
        point_clouds_ground[poly_i].push_back(point[2]);
        min_ground_elevation = std::min(min_ground_elevation, point[2]);
      } else if (point_class == 6) {
        if (GridTest(poly_grids[poly_i], pipoint)) {
          point_clouds.get<PointCollection&>(poly_i).push_back(point);
        }
      }
    }
    delete pipoint;
  }

  void compute_ground_elevation(float& ground_percentile) {
    // Compute median ground height per grid cell and store it for each polygon
    std::cout <<"Computing the median ground elevation per polygon..." << std::endl;
    for (auto& z_vec : point_clouds_ground) {
      // Median for odd-number of elements
      float ground_ele = min_ground_elevation;
      if (z_vec.size()!=0) {
        ground_ele = compute_percentile(z_vec, ground_percentile);
      } else {
        std::cout << "no ground pts found for polygon\n";
      }
      // Assign the median ground elevation to each polygon
      ground_elevations.push_back(ground_ele);
    }
  }
};


void LASInPolygonsNode::process() {
  auto& polygons = vector_input("polygons");

  auto& point_clouds = vector_output("point_clouds");
  auto& ground_elevations = vector_output("ground_elevations");

  PointsInPolygonsCollector pip_collector{polygons, point_clouds, ground_elevations, cellsize, buffer};

  LASreadOpener lasreadopener;
  lasreadopener.set_file_name(manager.substitute_globals(filepath).c_str());
  LASreader* lasreader = lasreadopener.open();
  
  if (!lasreader)
    return;

  while (lasreader->read_point()) {
    pip_collector.add_point(
      {
      float(lasreader->point.get_x()-(*manager.data_offset)[0]), 
      float(lasreader->point.get_y()-(*manager.data_offset)[1]),
      float(lasreader->point.get_z()-(*manager.data_offset)[2]) 
      }, 
      lasreader->point.get_classification()
    );
  }
  lasreader->close();
  delete lasreader;
  pip_collector.compute_ground_elevation(ground_percentile);
}

void EptInPolygonsNode::process()
{
  // Prepare inputs
  auto& polygons = vector_input("polygons");

  // Prepare outputs
  auto& point_clouds = vector_output("point_clouds");
  auto& ground_elevations = vector_output("ground_elevations");

  PointsInPolygonsCollector pip_collector{polygons, point_clouds, ground_elevations, cellsize, buffer};

  std::string ept_path = "ept://" + manager.substitute_globals(dirpath);
  {
    // This scope is for debugging
    pdal::EptReader reader;
    pdal::Options   options;
    options.add("filename", ept_path);
    reader.setOptions(options);
    // .preview() only reads the metadata
    const pdal::QuickInfo qi(reader.preview());
    std::cout << std::endl << "EPT Bounds:\t" << qi.m_bounds << std::endl;
    std::cout << "EPT Point Count:\t" << qi.m_pointCount << std::endl;
    std::cout << "EPT WKT:\t" << qi.m_srs.getWKT() << std::endl;
  }

  // TODO: The reader and filter init can go to a function that returns a
  //  const PointViewSet.
  auto pmin = pip_collector.completearea_bb.min();
  auto pmax = pip_collector.completearea_bb.max();
  pdal::BOX2D poly_bbox(
    pmin[0] - buffer + (*manager.data_offset)[0],
    pmin[1] - buffer + (*manager.data_offset)[1],
    pmax[0] + buffer + (*manager.data_offset)[0],
    pmax[1] + buffer + (*manager.data_offset)[1]
  );
  pdal::EptReader reader;
  {
    pdal::Options options;
    options.add("filename", ept_path);
    options.add("bounds", poly_bbox);
//      options.add("polygon", wkt.str());
    reader.setOptions(options);
  }
  pdal::RangeFilter range;
  {
    pdal::Options options;
    if (filter_limits.length() == 0) {
      std::cout << "Uh oh, PDAL Range filter cannot be empty. GOING TO CRASH!"
                << std::endl;
    }
    options.add("limits", filter_limits);
    range.setOptions(options);
    range.setInput(reader);
  }
  pdal::PointTable ept_table;
  range.prepare(ept_table);
  const pdal::PointViewSet pw_set(range.execute(ept_table));

  std::cout << "...PDAL Pipeline executed\n Continuing with point-in-polygon tests\n";

  for (const pdal::PointViewPtr& view : pw_set) {

    if (view->size() == 0) {
      std::cout << "0 points were found" << std::endl;
    }
    for (pdal::point_count_t p(0); p < view->size(); ++p) {
      float px = view->getFieldAs<float>(pdal::Dimension::Id::X, p) - (*manager.data_offset)[0];
      float py = view->getFieldAs<float>(pdal::Dimension::Id::Y, p) - (*manager.data_offset)[1];
      float pz = view->getFieldAs<float>(pdal::Dimension::Id::Z, p) - (*manager.data_offset)[2];
      int pclass = view->getFieldAs<int>(pdal::Dimension::Id::Classification, p);
      
      pip_collector.add_point({px, py, pz}, pclass);
    }
  }

  pip_collector.compute_ground_elevation(ground_percentile);
}

}