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
    vals_ = std::make_unique<std::vector<double>>();
    vals_->resize(dimx_*dimy_);
    counts_ = std::make_unique<std::vector<int16_t>>();
    counts_->resize(dimx_*dimy_);
  }
  Raster::Raster(const Raster& r)
  {
    cellSize_ = r.cellSize_;
    maxx_ = r.maxx_;
    minx_ = r.minx_;
    maxy_ = r.maxy_;
    miny_ = r.miny_;
    noDataVal_ = r.noDataVal_;
    dimx_ = (maxx_-minx_)/cellSize_ + 1;
    dimy_ = (maxy_-miny_)/cellSize_ + 1;
    vals_ = std::make_unique<std::vector<double>>(*r.vals_);
    counts_ = std::make_unique<std::vector<int16_t>>(*r.counts_);
  }
  void Raster::operator=(const Raster& r)
  {
    cellSize_ = r.cellSize_;
    maxx_ = r.maxx_;
    minx_ = r.minx_;
    maxy_ = r.maxy_;
    miny_ = r.miny_;
    noDataVal_ = r.noDataVal_;
    dimx_ = (maxx_-minx_)/cellSize_ + 1;
    dimy_ = (maxy_-miny_)/cellSize_ + 1;
    vals_ = std::make_unique<std::vector<double>>(*r.vals_);
    counts_ = std::make_unique<std::vector<int16_t>>(*r.counts_);
  }

  void Raster::prefill_arrays(alg a){
  if (a==MIN)
    noDataVal_ = 99999;
  else if (a==CNT)
    noDataVal_ = 0;
  else
    noDataVal_ = -99999;
    
  std::fill(vals_->begin(), vals_->end(), noDataVal_);
  std::fill(counts_->begin(), counts_->end(), 0);
  }

  bool Raster::add_point(double x, double y, double z, alg a)
  {
    bool first = (*vals_)[getLinearCoord(x,y)]==noDataVal_;
    if (a==MIN) {
      min(x,y,z);
    } else if (a==MAX) {
      max(x,y,z);
    } else if (a==AVG) {
      avg(x,y,z);
    } else if (a==CNT) {
      cnt(x,y);
    }
    return first;
  }

  inline void Raster::avg(double &x, double &y, double &val)
  {
    size_t c = getLinearCoord(x,y);
    (*vals_)[c]= ((*vals_)[c]*(*counts_)[c]+val)/((*counts_)[c]+1);
    ++(*counts_)[c];
  }

  inline void Raster::min(double &x, double &y, double &val)
  {
    size_t c = getLinearCoord(x,y);
    if ((*vals_)[c]>val) (*vals_)[c] = val;
  }

  inline void Raster::max(double &x, double &y, double &val)
  {
    size_t c = getLinearCoord(x,y);
    if ((*vals_)[c]<val) (*vals_)[c] = val;
  }

  inline void Raster::cnt(double &x, double &y)
  {
    size_t c = getLinearCoord(x,y);
    ++(*counts_)[c];
  }

  std::array<double,2> Raster::getColRowCoord(double x, double y) const
  {
    double r = (y-miny_) / cellSize_;
    double c = (x-minx_) / cellSize_;
    
    return {c,r};
  }

  size_t Raster::getRow(double x, double y) const
  {
    return static_cast<size_t>( floor((y-miny_) / cellSize_) );
  }
  size_t Raster::getCol(double x, double y) const
  {
    return static_cast<size_t>( floor((x-minx_) / cellSize_) );
  }

  size_t Raster::getLinearCoord(double x, double y) const
  {
    size_t r = static_cast<size_t>( floor((y-miny_) / cellSize_) );
    size_t c = static_cast<size_t>( floor((x-minx_) / cellSize_) );
    
    return r * dimx_ + c;
  }

  std::array<float,3> Raster::getPointFromRasterCoords(size_t col, size_t row) const
  {
    std::array<float,3> p;
    p[0] = minx_ + col*cellSize_ + cellSize_/2;
    p[1] = miny_ + row*cellSize_ + cellSize_/2;
    p[2] = (*vals_)[col+row*dimx_];
    return p;
  }

  double Raster::sample(double &x, double &y)
  {
    return (*vals_)[getLinearCoord(x,y)];
  }

  void Raster::set_val(size_t col, size_t row, double val) {
    (*vals_)[col+row*dimx_] = val;
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