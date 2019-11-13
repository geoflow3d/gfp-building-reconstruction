// Copyright (c) 2012, 2013
// Ravi Peters -- r.y.peters@tudelft.nl
// All rights reserved
// 
// This file is part of Surfonoi.
// 
// Surfonoi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Surfonoi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Surfonoi. If not, see <http://www.gnu.org/licenses/>.

#include "Raster.h"

namespace RasterTools {

  Raster::Raster(double cellsize, double min_x, double max_x, double min_y, double max_y):
      cellSize_(cellsize), minx_(min_x), maxx_(max_x), miny_(min_y), maxy_(max_y)
  {  
    dimx_ = (maxx_-minx_)/cellSize_ + 1;
    dimy_ = (maxy_-miny_)/cellSize_ + 1;
    vals_ = new double[dimx_*dimy_];
    counts_ = new int16_t[dimx_*dimy_];
  }

  Raster::~Raster()
  {
    delete[] vals_;
    delete[] counts_;
  }

  void Raster::prefill_arrays(alg a){
  if (a==MIN)
    noDataVal_ = 99999;
  else if (a==CNT)
    noDataVal_ = 0;
  else
    noDataVal_ = -99999;
    
  std::fill(vals_+0, vals_+dimx_*dimy_, noDataVal_);
  std::fill(counts_+0, counts_+dimx_*dimy_, 0);
  }

  void Raster::add_point(double x, double y, double z, alg a)
  {
    if (a==MIN) {
      min(x,y,z);
    } else if (a==MAX) {
      max(x,y,z);
    } else if (a==AVG) {
      avg(x,y,z);
    } else if (a==CNT) {
      cnt(x,y);
    }
  }

  inline void Raster::avg(double &x, double &y, double &val)
  {
    size_t c = getLinearCoord(x,y);
    vals_[c]= (vals_[c]*counts_[c]+val)/(counts_[c]+1);
    ++counts_[c];
  }

  inline void Raster::min(double &x, double &y, double &val)
  {
    size_t c = getLinearCoord(x,y);
    if (vals_[c]>val) vals_[c] = val;
  }

  inline void Raster::max(double &x, double &y, double &val)
  {
    size_t c = getLinearCoord(x,y);
    if (vals_[c]<val) vals_[c] = val;
  }

  inline void Raster::cnt(double &x, double &y)
  {
    size_t c = getLinearCoord(x,y);
    ++counts_[c];
  }

  inline std::array<double,2> Raster::getColRowCoord(double x, double y)
  {
    double r = (y-miny_) / cellSize_;
    double c = (x-minx_) / cellSize_;
    
    return {c,r};
  }

  size_t Raster::getLinearCoord(double &x, double &y)
  {
    size_t r = static_cast<size_t>( floor((y-miny_) / cellSize_) );
    size_t c = static_cast<size_t>( floor((x-minx_) / cellSize_) );
    
    return r * dimx_ + c;
  }

  std::array<float,3> Raster::getPointFromRasterCoords(size_t col, size_t row) 
  {
    std::array<float,3> p;
    p[0] = minx_ + col*cellSize_ + cellSize_/2;
    p[1] = miny_ + row*cellSize_ + cellSize_/2;
    double v = vals_[col+row*dimx_];
    if (v==noDataVal_)
      p[2] = 0;
    else
      p[2] = v;
    return p;
  }

  double Raster::sample(double &x, double &y)
  {
    return vals_[getLinearCoord(x,y)];
  }

  std::vector<Raster::point3d> Raster::rasterise_polygon(std::vector<point3d>& polygon) {
    // code adapted from http://alienryderflex.com/polygon_fill/
    int n_nodes, pixelX, pixelY, i, j, swap ;
    int n_vertices = polygon.size();
    std::vector<point3d> result;

    // perhaps we can specialise these to the bounding box of the polygon
    int IMAGE_TOP = 0, IMAGE_BOT = dimy_, IMAGE_LEFT = 0, IMAGE_RIGHT=dimx_;

    // Loop through the rows of the image.
    for (pixelY=IMAGE_TOP; pixelY<IMAGE_BOT; pixelY++) {
      std::vector<int> intersect_x; // vector to hold the x-coordinates where the scanline intersects the polygon

      // Build a list of nodes.
      n_nodes=0; j=n_vertices-1;
      for (i=0; i<n_vertices; i++) {
        auto pi = getColRowCoord((double)polygon[i][0], (double)polygon[i][1]);
        auto pj = getColRowCoord((double)polygon[j][0], (double)polygon[j][1]);
        // std::cerr << pi[0] << " " << pi[1] << "\n";
        // std::cerr << pj[0] << " " << pj[1] << "\n";
        if ( (pi[1]<(double) pixelY && pj[1]>=(double) pixelY)
        || (pj[1]<(double) pixelY && pi[1]>=(double) pixelY)) {
          intersect_x.push_back((int) (pi[0]+(pixelY-pi[1])/(pj[1]-pi[1])
          *(pj[0]-pi[0])));
          ++n_nodes;
        }
        j=i; 
      }

      // Sort the nodes, via a simple “Bubble” sort.
      i=0;
      while (i<n_nodes-1) {
        if (intersect_x[i]>intersect_x[i+1]) {
          swap=intersect_x[i]; intersect_x[i]=intersect_x[i+1]; intersect_x[i+1]=swap; if (i) i--; 
        } else {
          i++; 
        }
      }

      // Fill the pixels between node pairs.
      for (i=0; i<n_nodes; i+=2) {
        if  (intersect_x[i ]>=IMAGE_RIGHT) 
          break;
        if  (intersect_x[i+1]> IMAGE_LEFT ) {
          if (intersect_x[i ]< IMAGE_LEFT ) 
            intersect_x[i ]=IMAGE_LEFT ;
          if (intersect_x[i+1]> IMAGE_RIGHT) 
            intersect_x[i+1]=IMAGE_RIGHT;
          for (pixelX=intersect_x[i]; pixelX<intersect_x[i+1]; pixelX++) 
            result.push_back(getPointFromRasterCoords(pixelX,pixelY)); 
        }
      }
    }

    return result;
  }

  // void Raster::write(const char* WKGCS, alg a, void * dataPtr, const char* outFile)
  // {
  //   if( EQUALN(WKGCS, "EPSG:",5) ) {
  //     oSRS.importFromEPSG( atoi(WKGCS+5) );
  //   } else if (EQUALN(WKGCS, "EPSGA:",6)) {
  //     oSRS.importFromEPSGA( atoi(WKGCS+6) );
  //   }
  //   GDALAllRegister();
  //   GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
  //   GDALDataset *poDstDS;
  //   GDALDataType dataType;

  //   if (a == CNT)
  //     dataType = GDT_UInt16;
  //   else
  //     dataType = GDT_Float64;
    
  //   char **papszOptions = NULL;
  //   poDstDS = poDriver->Create( outFile, dimx_, dimy_, 1, dataType,
  //                papszOptions );
  //   double adfGeoTransform[6] = { minx_, cellSize_, 0, miny_, 0, cellSize_ };
  //   GDALRasterBand *poBand;
    
  //   poDstDS->SetGeoTransform( adfGeoTransform );
    
  //   //  std::cout << oSRS.SetWellKnownGeogCS( WKGCS );
  //   //  std::cout << pszSRS_WKT <<std::endl;
    
  //   char *pszSRS_WKT = NULL;
  //   oSRS.exportToWkt( &pszSRS_WKT );
  //   poDstDS->SetProjection( pszSRS_WKT );
  //   CPLFree( pszSRS_WKT );
    
  //   poBand = poDstDS->GetRasterBand(1);
  //   poBand->RasterIO( GF_Write, 0, 0, dimx_, dimy_,
  //           dataPtr, dimx_, dimy_, dataType, 0, 0 );
  //   poBand->SetNoDataValue(noDataVal);
  //   /* Once we're done, close properly the dataset */
  //   GDALClose( (GDALDatasetH) poDstDS );
  // }

}