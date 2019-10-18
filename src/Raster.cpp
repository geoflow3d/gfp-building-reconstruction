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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Surfonoi.  If not, see <http://www.gnu.org/licenses/>.

#include "Raster.h"

namespace RasterTools {

    Raster::Raster(double cellsize, double min_x, double max_x, double min_y, double max_y):
            cellSize(cellsize), minx(min_x), maxx(max_x), miny(min_y), maxy(max_y)
    {    
        dimx = (maxx-minx)/cellSize + 1;
        dimy = (maxy-miny)/cellSize + 1;
        vals = new double[dimx*dimy];
        counts = new int16_t[dimx*dimy];
    }

    Raster::~Raster()
    {
        delete[] vals;
        delete[] counts;
    }

    void Raster::prefill_arrays(alg a){
    if (a==MIN)
        noDataVal = 99999;
    else if (a==CNT)
        noDataVal = 0;
    else
        noDataVal = -99999;
        
    std::fill(vals+0, vals+dimx*dimy, noDataVal);
    std::fill(counts+0, counts+dimx*dimy, 0);
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
        vals[c]= (vals[c]*counts[c]+val)/(counts[c]+1);
        ++counts[c];
    }

    inline void Raster::min(double &x, double &y, double &val)
    {
        size_t c = getLinearCoord(x,y);
        if (vals[c]>val) vals[c] = val;
    }

    inline void Raster::max(double &x, double &y, double &val)
    {
        size_t c = getLinearCoord(x,y);
        if (vals[c]<val) vals[c] = val;
    }

    inline void Raster::cnt(double &x, double &y)
    {
        size_t c = getLinearCoord(x,y);
        ++counts[c];
    }

    size_t Raster::getLinearCoord(double &x, double &y)
    {
        size_t r = static_cast<size_t>( floor((y-miny) / cellSize) );
        size_t c = static_cast<size_t>( floor((x-minx) / cellSize) );
        
        return r * dimx + c;
    }

    std::array<float,3> Raster::getPointFromRasterCoords(size_t col, size_t row) 
    {
        std::array<float,3> p;
        p[0] = minx + col*cellSize + cellSize/2;
        p[1] = miny + row*cellSize + cellSize/2;
        double v = vals[col+row*dimx];
        if (v==noDataVal)
            p[2] = 0;
        else
            p[2] = v;
        return p;
    }

    double Raster::sample(double &x, double &y)
    {
        return vals[getLinearCoord(x,y)];
    }

    // void Raster::write(const char* WKGCS, alg a, void * dataPtr, const char* outFile)
    // {
    //     if( EQUALN(WKGCS, "EPSG:",5) ) {
    //         oSRS.importFromEPSG( atoi(WKGCS+5) );
    //     } else if (EQUALN(WKGCS, "EPSGA:",6)) {
    //         oSRS.importFromEPSGA( atoi(WKGCS+6) );
    //     }
    //     GDALAllRegister();
    //     GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    //     GDALDataset *poDstDS;
    //     GDALDataType dataType;

    //     if (a == CNT)
    //         dataType = GDT_UInt16;
    //     else
    //         dataType = GDT_Float64;
        
    //     char **papszOptions = NULL;
    //     poDstDS = poDriver->Create( outFile, dimx, dimy, 1, dataType,
    //                                papszOptions );
    //     double adfGeoTransform[6] = { minx, cellSize, 0, miny, 0, cellSize };
    //     GDALRasterBand *poBand;
        
    //     poDstDS->SetGeoTransform( adfGeoTransform );
        
    //     //    std::cout << oSRS.SetWellKnownGeogCS( WKGCS );
    //     //    std::cout << pszSRS_WKT <<std::endl;
        
    //     char *pszSRS_WKT = NULL;
    //     oSRS.exportToWkt( &pszSRS_WKT );
    //     poDstDS->SetProjection( pszSRS_WKT );
    //     CPLFree( pszSRS_WKT );
        
    //     poBand = poDstDS->GetRasterBand(1);
    //     poBand->RasterIO( GF_Write, 0, 0, dimx, dimy,
    //                      dataPtr, dimx, dimy, dataType, 0, 0 );
    //     poBand->SetNoDataValue(noDataVal);
    //     /* Once we're done, close properly the dataset */
    //     GDALClose( (GDALDatasetH) poDstDS );
    // }

}