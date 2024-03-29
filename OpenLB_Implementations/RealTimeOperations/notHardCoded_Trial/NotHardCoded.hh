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

 /* Important Note:
 *  This is just a trial for implementing OpenLB as a generalized program.
 *  Not fully functional!
 */

#ifndef NOT_HARD_CODED_HH
#define NOT_HARD_CODED_HH

#include "NotHardCoded.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using T = FLOATING_POINT_TYPE;

// Usage of the typedef enums must be added!
typedef enum {D3Q7, D3Q19, D3Q27} DescriptorType;
typedef enum {forced, nonForced} FlowType;
typedef enum {heated, nonHeated} HeatType;
typedef enum {bounceBack, local, interpolated, bouzidi, freeSlip, partialSlip} BoundaryType;
typedef enum {Cuboid3D, Cylinder3D, Circle3D} GeometryType;

// Parameters for the simulation setup 
DescriptorType descriptorType;
FlowType flowType;
HeatType heatType;
BoundaryType boundaryType;
GeometryType geometryType;

std::vector<std::string> pathSTL;
std::string filepath = "";
std::vector<T> converterProperties;
int tracker = 0;

template<typename T>
ModelMSO<T>::ModelMSO(std::vector<std::string> stringToStart)
{
    for(const std::string& str : stringToStart)
    {
        float floatValue;
        std::istringstream iss(str);
        if(iss >> floatValue)
            _converterProperties.push_back(floatValue);
    }

    _numberOfObjects = 0;
}

template<typename T>
ModelMSO<T>::~ModelMSO() {
}

template<typename T>
bool ModelMSO<T>::getGeometryData()
{
    bool RedisConnect = false;
    redisContext *redis = redisConnect("127.0.0.1", 6379);

    if (redis == NULL || redis->err) {
        std::cout << "Failed to connect to Redis: " << (redis ? redis->errstr : "Connection error") << std::endl;
        if (redis)
            redisFree(redis);
        return false;
    }

    redisReply *reply = (redisReply *)redisCommand(redis, "LRANGE geometryValues 0 -1");

    std::cout << "Number of elements: " << reply->elements << std::endl;

    if (reply != NULL && reply->type == REDIS_REPLY_ARRAY && reply->elements > 7) 
    {
        for (size_t i = 0; i < reply->elements -1; ++i) {
            _geometryString.push_back(reply->element[i]->str);
        }

        freeReplyObject(reply);
        RedisConnect = true;
    } 
    else 
    {
        std::cout << "Failed to retrieve strings from Redis." << std::endl;
    }

    redisFree(redis);

    return RedisConnect;
}

template<typename T>
bool ModelMSO<T>::getAndSetBoundaryData()
{
    bool RedisConnect = false;
    redisContext *redis = redisConnect("127.0.0.1", 6379);

    if (redis == NULL || redis->err) {
        std::cout << "Failed to connect to Redis: " << (redis ? redis->errstr : "Connection error") << std::endl;
        if (redis)
            redisFree(redis);
        return false;
    }

    redisReply *reply = (redisReply *)redisCommand(redis, "LRANGE boundaryValues 0 -1");

    std::cout << "Number of elements: " << reply->elements << std::endl;

    if (reply != NULL && reply->type == REDIS_REPLY_ARRAY && reply->elements > 7) 
    {
        for (size_t i = 0; i < reply->elements -1; ++i) {
            _boundaryProperties.push_back(std::stof(reply->element[i]->str));
        }

        freeReplyObject(reply);
        RedisConnect = true;
    } 
    else 
    {
        std::cout << "Failed to retrieve strings from Redis." << std::endl;
    }

    redisFree(redis);

    return RedisConnect;
}

template<typename T>
void ModelMSO<T>::setGeometryData()
{
    std::cout << "Geometry String: " << _geometryString << std::endl;
    float objectType = 0;
    for(int i = 0; i < _geometryString.size(); i++)
    {
        if(i == static_cast<int>(objectType))
        {
            if(_geometryString[i] == "Cube")
                _geometryProperties.push_back(0);
            if(_geometryString[i] == "Cylinder")
                _geometryProperties.push_back(1);
            if(_geometryString[i] == "Sphere")
                _geometryProperties.push_back(2);

            objectType += 8;
            _numberOfObjects++;
        }
        else
        {
            float floatValue = std::stof(_geometryString[i]);
            _geometryProperties.push_back(floatValue);
        }     
    }
}

template<typename T>
std::vector<T> ModelMSO<T>::getGeometryProperties()
{
    return _geometryProperties;
}

template<typename T>
void ModelMSO<T>::getIndicators(std::vector<std::shared_ptr<IndicatorF3D<T>>> indicators)
{
    _indicators = indicators;
}

template<typename T>
void ModelMSO<T>::prepareGeometry(SuperGeometry<T,3>& superGeometry,
                    std::vector<std::shared_ptr<IndicatorF3D<T>>> indicators)
{
    OstreamManager clout( std::cout,"prepareGeometry" );
    clout << "Prepare Geometry ..." << std::endl;

    superGeometry.rename( 0,2 );
    superGeometry.rename( 2,1,{1,1,1} );
    superGeometry.clean();

    _materialNumbers.push_back(1);
    _nonAdiabatic.push_back(1);

    for(int i = 0; i < _boundaryProperties.size(); i += 6)
    {
        _materialNumbers.push_back(static_cast<int>(_boundaryProperties[i]));

        if(static_cast<int>(_boundaryProperties[i+1]) != 2 )
        {
            _nonAdiabatic.push_back(static_cast<int>(_boundaryProperties[i]));
        }
    }

    std::cout << "Material Numbers: " << _materialNumbers << std::endl;

    // Set material number for cylinder

    std::cout << "Geomety Properties: " << _geometryProperties << std::endl;
    int j = 8;
    for(int i = 1; i < _geometryProperties.size() / 8; i++)
    {
        geometryType = static_cast<GeometryType>(_geometryProperties[j]);

        if(static_cast<int>(geometryType) == 0)
        {
            Vector<T, 3> origin3D(_geometryProperties[j+1], _geometryProperties[j+2], _geometryProperties[j+3]);
            Vector<T, 3> extend3D(_geometryProperties[j+4], _geometryProperties[j+5], _geometryProperties[j+6]);

            std::shared_ptr<IndicatorF3D<T>> cuboid = std::make_shared<IndicatorCuboid3D<T>>( extend3D, origin3D );
            indicators.push_back(cuboid);
            std::cout << "Cuboid Created..." << std::endl;
            clout << "Prep Geo for: " << _materialNumbers[i]  << std::endl;
            if(_materialNumbers[i+1] != 3 || _materialNumbers[i+1] != 4)
                superGeometry.rename( 1,_materialNumbers[i+1],cuboid );
            else
                superGeometry.rename( 2,_materialNumbers[i+1],cuboid );
        }
        else if(static_cast<int>(geometryType) == 1)
        {
            Vector<T, 3> center1(_geometryProperties[j+1], _geometryProperties[j+2], _geometryProperties[j+3]);
            Vector<T, 3> center2(_geometryProperties[j+4], _geometryProperties[j+5], _geometryProperties[j+6]);

            std::shared_ptr<IndicatorF3D<T>> cylinder = std::make_shared<IndicatorCylinder3D<T>>( center1, center2, _geometryProperties[j+7] );
            indicators.push_back(cylinder);
            std::cout << "Cylinder Created..." << std::endl;
            clout << "Prep Geo for: " << _materialNumbers[i]  << std::endl;
            if(_materialNumbers[i+1] != 3 || _materialNumbers[i+1] != 4)
                superGeometry.rename( 1,_materialNumbers[i+1],cylinder );
            else
                superGeometry.rename( 2,_materialNumbers[i+1],cylinder );
        }
        else if(static_cast<int>(geometryType) == 2)
        {
            Vector<T, 3> center(_geometryProperties[j+1], _geometryProperties[j+2], _geometryProperties[j+3]);
            Vector<T, 3> normal(_geometryProperties[j+4], _geometryProperties[j+5], _geometryProperties[j+6]);

            std::shared_ptr<IndicatorF3D<T>> sphere = std::make_shared<IndicatorSphere3D<T>>( center, _geometryProperties[j+7] );
            indicators.push_back(sphere);
            std::cout << "Sphere Created..." << std::endl;
            clout << "Prep Geo for: " << _materialNumbers[i]  << std::endl;
            if(_materialNumbers[i+1] != 3 || _materialNumbers[i+1] != 4)
                superGeometry.rename( 1,_materialNumbers[i+1],sphere );
            else
                superGeometry.rename( 2,_materialNumbers[i+1],sphere );
        }
        else
        {
            std::cout << "[Error createYourGeo]" << std::endl;
        }
        j += 8;
    }

    // Removes all not needed boundary voxels outside the surface
    superGeometry.clean();
    superGeometry.innerClean();
    superGeometry.checkForErrors();

    superGeometry.print();


    clout << "Prepare Geometry ... OK" << std::endl;
}

template<typename T>
void ModelMSO<T>::prepareLatticeThermal( SuperLattice<T,NSDESCRIPTOR>& nsLattice,
                    SuperLattice<T,TDESCRIPTOR>& adLattice,
                    ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                    SuperGeometry<T,3>& superGeometry)
{
    using namespace olb;
    using namespace olb::descriptors;
    using namespace olb::graphics;
    OstreamManager clout( std::cout,"prepareLattice" );
    clout << "Prepare Lattice ..." << std::endl;
    
    T Tomega  = converter.getLatticeThermalRelaxationFrequency();
    T NSomega = converter.getLatticeRelaxationFrequency();

    AnalyticalConst3D<T,T> T_cold(0);
    AnalyticalConst3D<T,T> T_hot(1);

    AnalyticalConst3D<T,T> rhoF( 1 );
    std::vector<T> velocityzero( 3,T( 0 ) );
    AnalyticalConst3D<T,T> uF( velocityzero );

    T smagoConst = 0.15;

    // Material=1 -->bulk dynamics
    auto bulkIndicator = superGeometry.getMaterialIndicator({1});
    nsLattice.template defineDynamics<ExternalTauEffLESForcedBGKdynamics<T,NSDESCRIPTOR>>(superGeometry, 1);
    adLattice.template defineDynamics<ExternalTauEffLESBGKadvectionDiffusionDynamics>(superGeometry.getMaterialIndicator({1,2,3,4,5,6,7,8,9,10}));
    nsLattice.template setParameter<collision::LES::Smagorinsky>(smagoConst);
    adLattice.template setParameter<collision::LES::Smagorinsky>(smagoConst);

    std::vector<T> poiseuilleForce( 3,T() );
    poiseuilleForce[1] = -9.81/converter.getConversionFactorVelocity()*converter.getConversionFactorTime();
    AnalyticalConst3D<T,T> force( poiseuilleForce );

    // Initialize force
    nsLattice.template defineField<FORCE>(superGeometry, 1, force);
    nsLattice.template defineField<FORCE>(superGeometry, 2, force);

    // Material --> Boundary Type
    /*
        *Add local boundary option
        *Heat Flux
    */
    for(int i = 0; i < _boundaryProperties.size(); i += 6)
    {
        if(static_cast<int>(_boundaryProperties[i+2]) != 2 )
        {
            if(static_cast<int>(_boundaryProperties[i+1]) == 0 )
            {
                setBounceBackBoundary(nsLattice, superGeometry, static_cast<int>(_boundaryProperties[i]));
            }
            else
            {
                //Change to Bouzidi
                setBounceBackBoundary(nsLattice, superGeometry, static_cast<int>(_boundaryProperties[i]));
            }
        }
        else
        {
            if(_boundaryProperties[i+3] == 1)
            {
                setInterpolatedVelocityBoundary(nsLattice, NSomega, superGeometry, _boundaryProperties[i]);
                std::vector<T> poiseuilleForce( 3,T() );
                poiseuilleForce[2] = converter.getCharLatticeVelocity();
                AnalyticalConst3D<T,T> force( poiseuilleForce );
                AnalyticalConst3D<T,T> rhoFatm( 1. );
                nsLattice.defineU( superGeometry, _boundaryProperties[i], force );
                nsLattice.iniEquilibrium(superGeometry, _boundaryProperties[i], rhoFatm, force);
                clout << "Velocity Boundary for: " << _boundaryProperties[i]  << std::endl;
            }

            if(_boundaryProperties[i+3] == 3)
            {
                setInterpolatedPressureBoundary(nsLattice, NSomega, superGeometry, _boundaryProperties[i]);
                AnalyticalConst3D<T,T> rhoFatm( 1. );
                nsLattice.defineRho( superGeometry, _boundaryProperties[i], rhoFatm );
                clout << "Pressure Boundary for: " << _boundaryProperties[i]  << std::endl;
            }
        }
    }

    nsLattice.template defineField<FORCE>(superGeometry, 2, force);

    nsLattice.defineRho( superGeometry, 1, rhoF);
    nsLattice.defineU( superGeometry, 1, uF );
    nsLattice.iniEquilibrium( superGeometry, 1, rhoF, uF );

    adLattice.defineRho(superGeometry, 1, T_cold);
    adLattice.iniEquilibrium(superGeometry, 1, T_cold, uF);

    // Material --> Advection
    for(int i = 0; i < _boundaryProperties.size(); i += 6)
    {
        if(static_cast<int>(_boundaryProperties[i+1]) != 2 )
        {
            if(static_cast<int>(_boundaryProperties[i+1]) == 0 )
            {
                setAdvectionDiffusionTemperatureBoundary(adLattice, Tomega, superGeometry, static_cast<int>(_boundaryProperties[i]));
                AnalyticalConst3D<T,T> T_Source(converter.getLatticeTemperature(_boundaryProperties[i+4]));
                adLattice.defineRho(superGeometry, static_cast<int>(_boundaryProperties[i]), T_Source);
                adLattice.iniEquilibrium(superGeometry, static_cast<int>(_boundaryProperties[i]), T_Source, uF);
                clout << "Temp Boundary for: " << _boundaryProperties[i]  << std::endl;
            }
            else
            {
                //Heat Flux Later
                setAdvectionDiffusionTemperatureBoundary(adLattice, Tomega, superGeometry, static_cast<int>(_boundaryProperties[i]));
                AnalyticalConst3D<T,T> T_Source(converter.getLatticeTemperature(_boundaryProperties[i+4]));
                adLattice.defineRho(superGeometry, static_cast<int>(_boundaryProperties[i]), T_Source);
                adLattice.iniEquilibrium(superGeometry, static_cast<int>(_boundaryProperties[i]), T_Source, uF);
            }
        }
    }

    AnalyticalConst3D<T,T> tauNS(1./NSomega);
    AnalyticalConst3D<T,T> tauAD(1./Tomega);

    nsLattice.template defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1,2,3,4,5,6,7,8,9,10}), tauNS );
    adLattice.template defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1,2,3,4,5,6,7,8,9,10}), tauAD );

    nsLattice.template setParameter<descriptors::OMEGA>(NSomega);
    adLattice.template setParameter<descriptors::OMEGA>(Tomega);

    /// Make the lattice ready for simulation
    clout << "initialize ..." << std::endl;
    nsLattice.initialize();
    adLattice.initialize();

    clout << "Prepare Lattice ... OK" << std::endl;
}

template<typename T>
void ModelMSO<T>::getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );

  //SuperVTMwriter3D<T> vtmWriter("notHardCoded");
  const int vtkIter  = converter.getLatticeTime( .5 );
  const int statIter = converter.getLatticeTime( .1 );

  if ( iT%statIter == 0 ) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    NSlattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

    if ( iT%vtkIter == 0 && iT > converter.getLatticeTime(0.) ) 
    {
        SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( NSlattice, converter );
        AnalyticalFfromSuperF3D<T> velocityField( velocity, true );
        SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);

        T latticeR[3] { };
      
        int j = 0;

        for (latticeR[2]=_geometryProperties[3]; latticeR[2] <= _geometryProperties[6]; latticeR[2] += 2.*_converterProperties[1]/_converterProperties[0]) 
        {
            for (latticeR[1]=_geometryProperties[2]; latticeR[1] <= _geometryProperties[5]; latticeR[1] += 2.*_converterProperties[1]/_converterProperties[0]) 
            {
                for (latticeR[0]=_geometryProperties[1]; latticeR[0] <= _geometryProperties[4]; latticeR[0] += 2.*_converterProperties[1]/_converterProperties[0]) 
                {
                    T vel[3] { };
                    velocityField(vel, latticeR);
                    
                    if(tracker ==  0)
                        visualVector.push_back(0.0f);
                    else
                        visualVector[j] = std::sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
                
                    j++;
                }
            }
        }

        redisContext *redis = redisConnect("127.0.0.1", 6379);

        // Convert float vector to a string
        std::stringstream ss;
        for (const auto &element : visualVector) {
            ss << std::fixed << std::setprecision(4) << element << ",";
        }
        std::string vectorString = ss.str();
        vectorString.pop_back(); // Remove the trailing comma

        // Set the float vector in Redis
        redisReply *reply = (redisReply *)redisCommand(redis, "SET visualVector %s", vectorString.c_str());
        if (reply == nullptr || redis->err) {
            std::cerr << "Failed to set float vector in Redis: " << redis->errstr << std::endl;
        }

        freeReplyObject(reply);

        // Disconnect from Redis
        redisFree(redis);

        tracker++;
    }
}

template<typename T>
void ModelMSO<T>::simulate()
{
    using namespace olb;
    using namespace olb::descriptors;
    using namespace olb::graphics;

    // converterProperties[0] = cL
    // converterProperties[1] = N
    // converterProperties[2] = latticeVel
    // converterProperties[3] = Ra
    // converterProperties[4] = Pr
    // converterProperties[5] = Viscosity
    // converterProperties[6] = Density
    // converterProperties[7] = 25.684e-3,
    // converterProperties[8] = 0.00341
    // converterProperties[9] = Tcold
    // converterProperties[10] = Thot

    if(converterProperties.length() > 7)
    {
        ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter
        (
            (T) _converterProperties[1]/_converterProperties[0], // physDeltaX
            (T) 0.1/(0.4)*_converterProperties[1]/_converterProperties[0], // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
            (T) _converterProperties[1],  // charPhysLength
            (T) .1,//13.28 15.126e-6 / cL * util::sqrt( Ra / Pr ), // charPhysVelocity
            (T) _converterProperties[2],
            (T) _converterProperties[3],
            (T) 25.684e-3,
            (T) _converterProperties[8] * 25.684e-3 / _converterProperties[2] / _converterProperties[3],
            (T) 0.00341,
            (T) _converterProperties[5],
            (T) _converterProperties[6]
        );

        std::cout << "Converter Created..." << std::endl;
        converter.print();

        Vector<T, 3> origin3D(_geometryProperties[1], _geometryProperties[2], _geometryProperties[3]);
        Vector<T, 3> extend3D(_geometryProperties[4], _geometryProperties[5], _geometryProperties[6]);
        IndicatorCuboid3D cuboid3D(extend3D, origin3D);

        std::cout << "Domain Cuboid Created..." << std::endl;

        const int noOfCuboids = 1;

        CuboidGeometry3D<T> cuboidGeometry( cuboid3D, converter.getPhysDeltaX(), noOfCuboids, "volume" );

        std::cout << "Cuboid Geometry Created..." << std::endl;

        HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

        std::cout << "Load Balancer Created..." << std::endl;

        SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

        std::cout << "Super Geometry Created..." << std::endl;

        prepareGeometry(superGeometry, _indicators);

        SuperLattice<T, NSDESCRIPTOR> NSlattice( superGeometry );
        SuperLattice<T, TDESCRIPTOR> ADlattice( superGeometry );

        prepareLatticeThermal( NSlattice, ADlattice, converter, superGeometry );

        const T maxPhysT = 100;
        T smagoConst = 0.15;

        T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() * converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

        const T preFactor = smagoConst*smagoConst
                        * descriptors::invCs2<T,NSDESCRIPTOR>()*descriptors::invCs2<T,NSDESCRIPTOR>()
                        * 2*util::sqrt(2);

        SuperLatticeCoupling coupling(
        SmagorinskyBoussinesqCoupling{},
        names::NavierStokes{}, NSlattice,
        names::Temperature{},  ADlattice);
        coupling.template setParameter<SmagorinskyBoussinesqCoupling::T0>(
        converter.getLatticeTemperature(_converterProperties[5]));
        coupling.template setParameter<SmagorinskyBoussinesqCoupling::FORCE_PREFACTOR>(
        boussinesqForcePrefactor * Vector<T,3>{0.0,1.0,0.0});
        coupling.template setParameter<SmagorinskyBoussinesqCoupling::SMAGORINSKY_PREFACTOR>(preFactor);
        coupling.template setParameter<SmagorinskyBoussinesqCoupling::PR_TURB>(0.87);
        coupling.template setParameter<SmagorinskyBoussinesqCoupling::OMEGA_NSE>(
        converter.getLatticeRelaxationFrequency());
        coupling.template setParameter<SmagorinskyBoussinesqCoupling::OMEGA_ADE>(
        converter.getLatticeThermalRelaxationFrequency());

        std::cout << "Starting simulation..." << std::endl;
        util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
        timer.start();

        for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) 
        {
        // === 6th Step: Collide and Stream Execution ===
        NSlattice.collideAndStream();
        ADlattice.collideAndStream();

        coupling.execute();

        // === 7th Step: Computation and Output of the Results ===
        getResults( converter, NSlattice, ADlattice, iT, superGeometry, timer );
        }

        timer.stop();
        timer.printSummary();
    }
    else
    {
        //Non Thermal
    }

    
}

#endif