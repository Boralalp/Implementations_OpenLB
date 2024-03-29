/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2006-2019 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */
 
#ifndef NOT_HARD_CODED_H
#define NOT_HARD_CODED_H

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

#include <string>
#include <iostream>
#include <ios>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <hiredis/hiredis.h>
#include <chrono>
#include <thread>

    
template<typename T>
class ModelMSO {
private:
    std::vector<T> _converterProperties;
    std::vector<std::string> _geometryString;
    std::vector<T> _boundaryProperties;
    std::vector<int> _materialNumbers;
    std::vector<int> _nonAdiabatic;
    std::vector<T> _geometryProperties;
    int _numberOfObjects;
    int noOfCuboids;
    std::vector<std::shared_ptr<IndicatorF3D<T>>> _indicators;
    std::vector<T> visualVector;

public:
    using NSDESCRIPTOR = D3Q19<FORCE,TAU_EFF>;
    using TDESCRIPTOR = D3Q7<VELOCITY,TAU_EFF,tag::MRT>;

    ModelMSO(std::vector<std::string> stringToStart);
    virtual ~ModelMSO();

    bool getGeometryData();
    bool getAndSetBoundaryData();
    void setGeometryData();

    std::vector<T> getGeometryProperties();

    void getIndicators(std::vector<std::shared_ptr<IndicatorF3D<T>>> indicators);

    void prepareGeometry(SuperGeometry<T, 3>& superGeometry, std::vector<std::shared_ptr<IndicatorF3D<T>>> indicators);

    void prepareLatticeThermal(SuperLattice<T, NSDESCRIPTOR>& nsLattice,
                            SuperLattice<T, TDESCRIPTOR>& adLattice,
                            ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>& converter,
                            SuperGeometry<T, 3>& superGeometry);

    void createIndicator();

    void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer );

    void simulate();
};
#endif