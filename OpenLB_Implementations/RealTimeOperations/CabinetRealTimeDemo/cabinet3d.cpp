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
 
#include "olb3D.h"
#include "olb3D.hh"
 
using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
 
using T = FLOATING_POINT_TYPE;
 
using NSDESCRIPTOR = D3Q27<FORCE,TAU_EFF>;
using TDESCRIPTOR = D3Q27<VELOCITY,TAU_EFF>;
 
#define BOUZIDI
 
#include <string>
#include <iostream>
#include <ios>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <hiredis/hiredis.h>
#include <chrono>
#include <thread>
 
int N = 60;              // resolution of the model
const T cL = .6;         // char Length
const T maxPhysT = 200;  // max. simulation time in s, SI unit
const T latticeVel = 0.1;
 
const T Thot = 320;      // max Temperature
const T Tcold = 300;     // min Temperature
 
const T smagoConst = 0.2;

T L = cL/N;
T maxVelocity = 0;

std::vector<float> _converterProperties;
std::vector<std::string> _geometryData;
std::vector<T> visualVectorVelocity;
std::vector<T> visualVectorTemperature;
std::vector<std::string> fileName;

std::vector<std::string> _oldSimulationData = {"60", "2.25", "0.04", "300", "300", "0", "0"};
std::vector<std::string> _changedSimulationData = {"60", "2.25", "0.04", "300", "300", "0", "0"};
bool changeInSimulation = false;
bool firstTime = true;
 
// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T,TDESCRIPTOR>& ADlattice,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                     SuperGeometry<T,3>& superGeometry,
                     AnalyticalFfromSuperF3D<T>& fieldU,
                     AnalyticalFfromSuperF3D<T>& fieldT)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
 
  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();
 
  clout << "Tomega: " << Tomega << std::endl;
  clout << "Tomega: " << NSomega << std::endl;
 
  AnalyticalConst3D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst3D<T,T> T_hot(converter.getLatticeTemperature(Thot));
 
  AnalyticalConst3D<T,T> T_hotO1(converter.getLatticeTemperature(315.));
  AnalyticalConst3D<T,T> T_hotO2(converter.getLatticeTemperature(320.0));

  AnalyticalConst3D<T,T>rhoF( 1 );
  std::vector<T> velocityzero( 3,T( 0 ) );
  AnalyticalConst3D<T,T>uF( velocityzero );
 
  using ForcedSpecialdynamics = dynamics::Tuple<T, TDESCRIPTOR, momenta::BulkTuple, equilibria::SecondOrder, collision::OmegaFromCellTauEff<collision::BGK>, AdvectionDiffusionExternalVelocityCollision >;
 
  NSlattice.defineDynamics<ExternalTauEffLESForcedBGKdynamics<T,NSDESCRIPTOR>>(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5, 6, 7}));
  ADlattice.defineDynamics<ForcedSpecialdynamics>(superGeometry.getMaterialIndicator({1, 2, 4, 5, 6}));
  NSlattice.setParameter<collision::LES::Smagorinsky>(smagoConst);
  ADlattice.setParameter<collision::LES::Smagorinsky>(smagoConst);
 
  std::vector<T> gravityForce( 3,T() );
  gravityForce[1] = -9.81/converter.getConversionFactorVelocity()*converter.getConversionFactorTime();
  AnalyticalConst3D<T,T> gravity3D( gravityForce );
 
  // Initialize gravitiy force
  NSlattice.defineField<FORCE>(superGeometry, 1, gravity3D);
  NSlattice.defineField<FORCE>(superGeometry, 2, gravity3D);
 
  setLocalVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, NSomega, superGeometry.getMaterialIndicator({2,4,5}));

  setInterpolatedVelocityBoundary(NSlattice, NSomega, superGeometry, 3);
  setInterpolatedVelocityBoundary(NSlattice, NSomega, superGeometry, 6);

  setAdvectionDiffusionTemperatureBoundary(ADlattice, Tomega, superGeometry, 2);
  setAdvectionDiffusionTemperatureBoundary(ADlattice, Tomega, superGeometry, 4);
  setAdvectionDiffusionTemperatureBoundary(ADlattice, Tomega, superGeometry, 6);
 
  // Initialize all values of distribution functions to their local equilibrium
  NSlattice.defineRho( superGeometry, 1, rhoF);
  NSlattice.defineU( superGeometry, 1, fieldU );
  NSlattice.iniEquilibrium( superGeometry, 1, rhoF, fieldU );
  if(firstTime)
  {
    ADlattice.defineRho(superGeometry, 1, T_cold);
    ADlattice.iniEquilibrium(superGeometry, 1, T_cold, uF);
    firstTime = false;
  }
  else
  {
    NSlattice.defineU( superGeometry, 1, fieldU );
    NSlattice.iniEquilibrium( superGeometry, 1, rhoF, fieldU );

    AnalyticalConst3D<T,T> TconvertHot(Thot);
    AnalyticalConst3D<T,T> TconvertCold(Tcold);

    ADlattice.defineRho(superGeometry, 1, ((fieldT-TconvertCold)/(TconvertHot-TconvertCold))*(T_hot-T_cold)+T_cold);
    ADlattice.iniEquilibrium(superGeometry, 1, ((fieldT-TconvertCold)/(TconvertHot-TconvertCold))*(T_hot-T_cold)+T_cold, fieldU);
  }
  ADlattice.defineRho(superGeometry, 2, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 2, T_cold, uF);

  ADlattice.defineRho(superGeometry, 4, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 4, T_cold, uF);
  ADlattice.defineRho(superGeometry, 5, T_hotO1);
  ADlattice.iniEquilibrium(superGeometry, 5, T_hotO1, uF);
  ADlattice.defineRho(superGeometry, 6, T_hotO2);
  ADlattice.iniEquilibrium(superGeometry, 6, T_hotO2, uF);
 
  AnalyticalConst3D<T,T> tauNS(1./NSomega);
  AnalyticalConst3D<T,T> tauAD(1./Tomega);
 
  NSlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3, 4, 5, 6, 7}), tauNS );
  ADlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 4, 5, 6}), tauAD );
 
  NSlattice.setParameter<descriptors::OMEGA>(NSomega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);
 
  /// Make the lattice ready for simulation
  clout << "initialize ..." << std::endl;
  NSlattice.initialize();
  ADlattice.initialize();
 
  clout << "Prepare Lattice ... OK" << std::endl;
}
 
// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                        SuperLattice<T, TDESCRIPTOR>& ADlattice,
                        ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                        int iT,
                        SuperGeometry<T,3>& superGeometry)
{
 
  OstreamManager clout( std::cout,"setBoundaryValues" );
 
  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( 0.4 );
  int iTupdate = 100;
  int updateAfterInitialStart = converter.getLatticeTime( .2 );
 
  if ( iT%iTupdate==0 && iT<= iTmaxStart )
  {
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);
 
    // Smooth start curve, polynomial
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );
 
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );

    maxVelocity = converter.getCharLatticeVelocity()/2 + converter.getCharLatticeVelocity()/2 *frac[0];//converter.getLatticeVelocity(0.01725)*frac[0];
 
    std::vector<T> fanForce_01( 3,T() );
    fanForce_01[0] = -maxVelocity;
    std::vector<T> fanForce_02( 3,T() );
    fanForce_02[1] = maxVelocity/1.52;
 
    AnalyticalConst3D<T,T> forceApplied_01( fanForce_01 );
    AnalyticalConst3D<T,T> forceApplied_02( fanForce_02 );

    AnalyticalConst3D<T,T> rhoFatm( 1. );
    AnalyticalConst3D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
 
    NSlattice.defineU( superGeometry, 3, forceApplied_01 );
    NSlattice.iniEquilibrium(superGeometry, 3, rhoFatm, forceApplied_01);
    NSlattice.defineU( superGeometry, 6, forceApplied_02 );
    NSlattice.iniEquilibrium(superGeometry, 6, rhoFatm, forceApplied_02);

    clout << "Max Velocity = " << maxVelocity << std::endl;
 
    NSlattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    ADlattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if ( iT%updateAfterInitialStart==0 && iT > iTmaxStart )
  {
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    maxVelocity = converter.getLatticeVelocity(std::stof(_changedSimulationData[1]));

    std::vector<T> fanForce_01( 3,T() );
    fanForce_01[0] = -maxVelocity;
    std::vector<T> fanForce_02( 3,T() );
    fanForce_02[1] = maxVelocity/1.52;
 
    AnalyticalConst3D<T,T> forceApplied_01( fanForce_01 );
    AnalyticalConst3D<T,T> forceApplied_02( fanForce_02 );

    AnalyticalConst3D<T,T> rhoFatm( 1. );
    AnalyticalConst3D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
 
    NSlattice.defineU( superGeometry, 3, forceApplied_01 );
    NSlattice.iniEquilibrium(superGeometry, 3, rhoFatm, forceApplied_01);
    NSlattice.defineU( superGeometry, 6, forceApplied_02 );
    NSlattice.iniEquilibrium(superGeometry, 6, rhoFatm, forceApplied_02);

    AnalyticalConst3D<T,T> T_hotO2(converter.getLatticeTemperature(std::stof(_changedSimulationData[3])));
    AnalyticalConst3D<T,T> T_hotO1(converter.getLatticeTemperature(std::stof(_changedSimulationData[4])));
    
    ADlattice.defineRho(superGeometry, 5, T_hotO1);
    ADlattice.defineRho(superGeometry, 6, T_hotO2);
 
    clout << "Max Velocity = " << maxVelocity << std::endl;
 
    NSlattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    ADlattice.setProcessingContext(ProcessingContext::Simulation);
  }
}

void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 IndicatorTranslate3D<T>& domainSTLreader )
{
  OstreamManager clout( std::cout,"getResults" );

  const int vtkIter  = converter.getLatticeTime( .2 );
  const int statIter = converter.getLatticeTime( .1 );

  Vector<T,3> startPoint;
  startPoint = domainSTLreader.getMin();
  Vector<T,3> endPoint;
  endPoint = domainSTLreader.getMax();

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
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( NSlattice, converter );
    AnalyticalFfromSuperF3D<T> velocityField( velocity, true );
    SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
    AnalyticalFfromSuperF3D<T> temperatureField( temperature, true );

    clout << "Vector will be filled..." << std::endl;

    T latticeR[3] { };
  
    visualVectorVelocity.clear();
    visualVectorTemperature.clear();

    for (latticeR[2]=0; latticeR[2] <= 0.599; latticeR[2] += std::stof(_changedSimulationData[2])) 
    {
        for (latticeR[1]=0; latticeR[1] <= 2.399; latticeR[1] += std::stof(_changedSimulationData[2])) 
        {
            for (latticeR[0]=0; latticeR[0] <= 0.599; latticeR[0] += std::stof(_changedSimulationData[2])) 
            {
                T vel[3] { };
                T temp[3] { };
                velocityField(vel, latticeR);
                temperatureField(temp, latticeR);
                
                visualVectorVelocity.push_back(std::sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]));
                visualVectorTemperature.push_back(std::sqrt(temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]));
            }
        }
    }

    clout << "Size: " << visualVectorVelocity.size() << std::endl;
    clout << "Size: " << visualVectorTemperature.size() << std::endl;

    clout << "Vector Filled... Sending it to Unity" << std::endl;

    redisContext *redis = redisConnect("127.0.0.1", 6379);

    // Convert float vector to a string
    std::stringstream ss;
    if(std::stoi(_changedSimulationData[5]) == 0)
    {
      for (const auto &element : visualVectorVelocity) {
          ss << std::fixed << std::setprecision(4) << element << ",";
      }
    }
    else
    { 
      for (const auto &element : visualVectorTemperature) {
          ss << std::fixed << std::setprecision(4) << element << ",";
      }
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
  }
}

void readConverterValues()
{
  redisContext *redis = redisConnect("127.0.0.1", 6379);

  redisReply *response = (redisReply *)redisCommand(redis, "FLUSHALL");
  freeReplyObject(response);
  redisFree(redis);

  bool RedisStartValues = true;


  while(RedisStartValues)
  {
    redisContext *redis = redisConnect("127.0.0.1", 6379);
    redisReply *reply = (redisReply *)redisCommand(redis, "LRANGE converterValues 0 -1");

    if (reply != NULL && reply->type == REDIS_REPLY_ARRAY && reply->elements > 7) 
    {
      _converterProperties.clear();
      for (size_t i = 0; i < reply->elements; ++i) {
          _converterProperties.push_back(std::stof(reply->element[i]->str));
      }

      if(_converterProperties.size() > 7)
      {
        RedisStartValues = false;

        std::cout << "Start Properties: " << _converterProperties << std::endl;
      }
      freeReplyObject(reply);
    } 
    else 
    {
        std::cout << "Failed to retrieve strings from Redis in readConverterValues." << std::endl;
    }

    redisFree(redis);
    std::this_thread::sleep_for(std::chrono::seconds(2));
  }
}

bool readGeometryData()
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
      _geometryData.clear();
      for (size_t i = 0; i < reply->elements -1; ++i) {
          _geometryData.push_back(reply->element[i]->str);
      }
      std::cout << "Geometry Properties: " << _geometryData << std::endl;
      freeReplyObject(reply);
      RedisConnect = true;
    } 
    else 
    {
        std::cout << "Failed to retrieve strings from Redis in readGeometryData." << std::endl;
    }

    redisFree(redis);

    return RedisConnect;
}

bool changeInSimulationData()
{
    bool RedisConnect = false;
    redisContext *redis = redisConnect("127.0.0.1", 6379);

    if (redis == NULL || redis->err) {
        std::cout << "Failed to connect to Redis: " << (redis ? redis->errstr : "Connection error") << std::endl;
        if (redis)
            redisFree(redis);
        return false;
    }

    redisReply *reply = (redisReply *)redisCommand(redis, "LRANGE newProperties 0 -1");

    std::cout << "Number of elements for change in the simulation: " << reply->elements << std::endl;

    if (reply != NULL && reply->type == REDIS_REPLY_ARRAY && reply->elements > 6) 
    {
      _oldSimulationData = _changedSimulationData;
      _changedSimulationData.clear();
      for (size_t i = 0; i < reply->elements; ++i) {
          _changedSimulationData.push_back(reply->element[i]->str);
      }
      std::cout << "Changed Simulation Properties: " << _changedSimulationData << std::endl;
      freeReplyObject(reply);

      if(std::abs(std::stof(_oldSimulationData[0]) - std::stof(_changedSimulationData[0])) > 0.01 || std::abs(std::stof(_oldSimulationData[1]) - std::stof(_changedSimulationData[1]))  > 0.01 || std::abs(std::stof(_oldSimulationData[6]) - std::stof(_changedSimulationData[6]))  > 0.01)
      {
        RedisConnect = true;
        std::cout << "Old Simulation Properties: " << _oldSimulationData << std::endl;
      }
    } 
    else 
    {
        std::cout << "Failed to retrieve strings from Redis in changeInSimulationData." << std::endl;
    }

    redisFree(redis);

    return RedisConnect;
}

void prepareGeometry( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                      SuperGeometry<T,3>& superGeometry,
                      IndicatorF3D<T>& indicator,
                      IndicatorTranslate3D<T>& domainSTLreader,
                      std::vector<std::string> filePaths )
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2,indicator);
  superGeometry.rename(2,1,domainSTLreader);

  clout << "Domain Generated..." << std::endl;

  for (int i = 0; i < filePaths.size(); i++) 
  {
    clout << "Reading File: " << filePaths[i]<< std::endl;
    STLreader<T> stlReader(filePaths[i], L, .001, 1, false );
    std::array<T,3> shift;
    Vector<T,3> shiftVector;
    shiftVector = stlReader.getMin();
    shift[0] = -shiftVector[0];
    shift[1] = -shiftVector[1];
    shift[2] = -shiftVector[2];
    IndicatorTranslate3D<T>stlReaderOrigin(shift,stlReader);

    //Object is at P(0,0,0)
    clout << "Object at P(0,0,0)" << std::endl;
    
    auto targetSTL = std::find(_geometryData.begin(), _geometryData.end(), fileName[i]);
    if (targetSTL != _geometryData.end()) 
    {
        // Element found, get its index
        int index = std::distance(_geometryData.begin(), targetSTL);

        clout << "Index of the Object is " << index << std::endl;

        T xPosition = std::stof(_geometryData[index+1]);
        T yPosition = std::stof(_geometryData[index+2]);
        T zPosition = std::stof(_geometryData[index+3]);

        shift[0] = xPosition;
        shift[1] = yPosition;
        shift[2] = zPosition;
        IndicatorTranslate3D<T> stlReaderShifted(shift,stlReaderOrigin);

        if(fileName[i] == "chassis.stl" || fileName[i] == "plate_all.stl")
          superGeometry.rename(1, 2, stlReaderShifted);
        else if(fileName[i] == "bottomLeftInlet.stl" || fileName[i] == "topLeftOutlet.stl" || fileName[i] == "topRightOutlet.stl")
          superGeometry.rename(1, 4, stlReaderShifted);
        else 
          superGeometry.rename(1, 5, stlReaderShifted);

        clout << "Rename opeartion on superGeometry Done..." << std::endl;
    }
  }

  clout << "Generating Fans..." << std::endl;
  std::string inletObject = "inletFan";
  std::string outletObject = "outletFan";

  auto targetSTLinlet = std::find(_geometryData.begin(), _geometryData.end(), inletObject);
  auto targetSTLoutlet = std::find(_geometryData.begin(), _geometryData.end(), outletObject);
  if (targetSTLinlet != _geometryData.end()) 
  {
      // Fan Inlet Generator
      int index = std::distance(_geometryData.begin(), targetSTLinlet);

      T xPosition = std::stof(_geometryData[index+1]);
      T yPosition = std::stof(_geometryData[index+2]);
      T zPosition = std::stof(_geometryData[index+3]);
      T xLength = std::stof(_geometryData[index+4]);
      T radius = std::stof(_geometryData[index+5]);
      
      Vector<T, 3> center0(xPosition - xLength, yPosition, zPosition);
      Vector<T, 3> center1(xPosition + xLength, yPosition, zPosition);
      IndicatorCylinder3D<T> pipe(center0, center1, radius);

      superGeometry.rename(2, 3, pipe);
  }

  clout << "Inlet Fan Done..." << std::endl;

  if (targetSTLoutlet != _geometryData.end()) 
  {
      // Fan Outlet Generator
      int index = std::distance(_geometryData.begin(), targetSTLoutlet);

      T xPosition = std::stof(_geometryData[index+1]);
      T yPosition = std::stof(_geometryData[index+2]);
      T zPosition = std::stof(_geometryData[index+3]);
      T xLength = std::stof(_geometryData[index+4]);
      T yLength = std::stof(_geometryData[index+5]);
      T zLength = std::stof(_geometryData[index+6]);
      
      Vector<T, 3> origin(xPosition, yPosition, zPosition);
      Vector<T, 3> extend(xLength, yLength, zLength);
      IndicatorCuboid3D<T> cuboid(extend, origin);

      superGeometry.rename(2, 6, cuboid);
  }

  clout << "Outlet Fan Done..." << std::endl;

  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.print();
}
 
int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
 
  // The Initial Parameters for the simulation setup
  readConverterValues();

  std::shared_ptr<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>> converter = std::make_shared<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>>(
      (T) _converterProperties[1]/_converterProperties[0], // physDeltaX
      (T) 0.1/(_converterProperties[8] * 1.75)*_converterProperties[1]/_converterProperties[0], // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
      (T) _converterProperties[1],  // charPhysLength
      (T) _converterProperties[8],  // charPhysVelocity
      (T) _converterProperties[2],
      (T) _converterProperties[3],
      (T) 25.684e-3,
      (T) _converterProperties[7] * 25.684e-3 / _converterProperties[2] / _converterProperties[3],
      (T) 0.00341,
      (T) _converterProperties[5],
      (T) _converterProperties[6]
  );

  // The Initial Parameters for the simulatio
  converter->print();
  converter->write("cabinet3d");

  N = _converterProperties[0];
  L = _converterProperties[1]/_converterProperties[0];

  _oldSimulationData[0] = std::to_string(N);
  _changedSimulationData[0] = std::to_string(N);

  _oldSimulationData[1] = std::to_string(_converterProperties[8]);
  _changedSimulationData[1] = std::to_string(_converterProperties[8]);

  bool getGeometryData = false;
  while(!getGeometryData)
  {
    getGeometryData = readGeometryData();
    std::this_thread::sleep_for(std::chrono::seconds(2));
  }

  std::string folder_name = "../stl";
  std::string file_extension = ".stl"; 
  std::string domainPath = folder_name+"/Domain.stl";
  std::string domainFile = "Domain.stl";
  std::vector<std::string> filePaths;
  
  for (const auto& entry : std::filesystem::directory_iterator(folder_name)) {
      if (std::filesystem::is_regular_file(entry) && entry.path().extension() == file_extension) {
          filePaths.push_back(entry.path());
          fileName.push_back(entry.path().filename());
          clout << "File Name: " << entry.path().filename() << std::endl;
          clout << "File Path: " << entry.path() << std::endl;
      }
  }

  auto target = std::find(filePaths.begin(), filePaths.end(), domainPath);
  filePaths.erase(target);
  auto targetFile = std::find(fileName.begin(), fileName.end(), domainFile);
  fileName.erase(targetFile);

  // Find Domain.stl index
  auto targetSTL = std::find(_geometryData.begin(), _geometryData.end(), domainFile);
  int indexDomain = 0;

  if (targetSTL != _geometryData.end()) 
  {
      // Element found, get its index
      indexDomain = std::distance(_geometryData.begin(), targetSTL);

      clout << "Index of the Domain is " << indexDomain << std::endl;
  }


  std::shared_ptr<STLreader<T>> stlReader = std::make_shared<STLreader<T>>( folder_name+"/Domain.stl", L, .001, 1, false  );

  // Move Domain to Origin
  std::array<T,3> shift;
  Vector<T,3> shiftVector;
  shiftVector = stlReader->getMin();
  shift[0] = -shiftVector[0];
  shift[1] = -shiftVector[1];
  shift[2] = -shiftVector[2];
  std::shared_ptr<IndicatorTranslate3D<T>> stlReaderOrigin = std::make_shared<IndicatorTranslate3D<T>>( shift,*stlReader );
 
  // Move Domain to Desired Position
  T xPosition = std::stof(_geometryData[indexDomain + 1]);
  T yPosition = std::stof(_geometryData[indexDomain + 2]);
  T zPosition = std::stof(_geometryData[indexDomain + 3]);
 
  shift[0] = xPosition;
  shift[1] = yPosition;
  shift[2] = zPosition;
  std::shared_ptr<IndicatorTranslate3D<T>> stlReaderShifted = std::make_shared<IndicatorTranslate3D<T>>( shift,*stlReaderOrigin );
  std::shared_ptr<IndicatorLayer3D<T>> extendedDomain = std::make_shared<IndicatorLayer3D<T>>( *stlReaderShifted, L );

  int noOfCuboids = 1;
  std::shared_ptr<CuboidGeometry3D<T>> cuboidGeometry = std::make_shared<CuboidGeometry3D<T>>( *extendedDomain, L, noOfCuboids, "volume" );
  std::shared_ptr<HeuristicLoadBalancer<T>> loadBalancer = std::make_shared<HeuristicLoadBalancer<T>>( *cuboidGeometry );
  std::shared_ptr<SuperGeometry<T,3>> superGeometry = std::make_shared<SuperGeometry<T,3>>( *cuboidGeometry, *loadBalancer );

  prepareGeometry( *converter, *superGeometry, *extendedDomain, *stlReaderShifted, filePaths );

  // === 3rd Step: Prepare Lattice ===
  std::shared_ptr<SuperLattice<T, NSDESCRIPTOR>> NSlattice = std::make_shared<SuperLattice<T, NSDESCRIPTOR>>(*superGeometry);
  std::shared_ptr<SuperLattice<T, TDESCRIPTOR>> ADlattice = std::make_shared<SuperLattice<T, TDESCRIPTOR>>(*superGeometry);

  SuperLatticeVelocity3D<T, NSDESCRIPTOR> velocity( *NSlattice );
  AnalyticalFfromSuperF3D<T> fieldVelocity(velocity, true);
  SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(*ADlattice, *converter);
  AnalyticalFfromSuperF3D<T> fieldTemperature(temperature, true);
    
  prepareLattice(  *NSlattice, *ADlattice, *converter, *superGeometry, fieldVelocity, fieldTemperature);

  T boussinesqForcePrefactor = 9.81 / converter->getConversionFactorVelocity() * converter->getConversionFactorTime() * converter->getCharPhysTemperatureDifference() * converter->getPhysThermalExpansionCoefficient();

  const T preFactor = smagoConst*smagoConst
                    * descriptors::invCs2<T,NSDESCRIPTOR>()*descriptors::invCs2<T,NSDESCRIPTOR>()
                    * 2*util::sqrt(2);
  label:

  SuperLatticeCoupling coupling(
    SmagorinskyBoussinesqCoupling{},
    names::NavierStokes{}, *NSlattice,
    names::Temperature{},  *ADlattice);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::T0>(
    converter->getLatticeTemperature(Tcold));
  coupling.setParameter<SmagorinskyBoussinesqCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,3>{0.0,1.0,0.0});
  coupling.setParameter<SmagorinskyBoussinesqCoupling::SMAGORINSKY_PREFACTOR>(preFactor);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::PR_TURB>(0.87);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_NSE>(
    converter->getLatticeRelaxationFrequency());
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_ADE>(
    converter->getLatticeThermalRelaxationFrequency());

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter->getLatticeTime( maxPhysT ), superGeometry->getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter->getLatticeTime( maxPhysT ); ++iT ) 
  {
    setBoundaryValues( *NSlattice, *ADlattice, *converter, iT, *superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    NSlattice->collideAndStream();
    ADlattice->collideAndStream();

    coupling.execute();

    // === 7th Step: Computation and Output of the Results ===
    getResults( *converter, *NSlattice, *ADlattice, iT, *superGeometry, timer, *stlReaderOrigin );

    const int vtkIter  = converter->getLatticeTime( .2 );

    if(iT%vtkIter == 0)
    {
      changeInSimulation = changeInSimulationData();
      
      if(changeInSimulation)
      {
        changeInSimulation = false;

        clout << "There was a change in the simulation. Re-generating lattices." << std::endl;

        readConverterValues();

        N = _converterProperties[0];
        L = _converterProperties[1]/_converterProperties[0];

        std::shared_ptr<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>> newConverter = std::make_shared<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>>(
            (T) _converterProperties[1]/_converterProperties[0], // physDeltaX
            (T) latticeVel/(_converterProperties[8] * 1.75)*_converterProperties[1]/_converterProperties[0], // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
            (T) _converterProperties[1],  // charPhysLength
            (T) _converterProperties[8],  // charPhysVelocity
            (T) _converterProperties[2],
            (T) _converterProperties[3],
            (T) 25.684e-3,
            (T) _converterProperties[7] * 25.684e-3 / _converterProperties[2] / _converterProperties[3],
            (T) 0.00341,
            (T) _converterProperties[5],
            (T) _converterProperties[6]
        );

        std::shared_ptr<STLreader<T>> newStlReader = std::make_shared<STLreader<T>>( folder_name+"/Domain.stl", L, .001, 1, false  );
        std::array<T,3> shift;
        Vector<T,3> shiftVector;
        shiftVector = stlReader->getMin();
        shift[0] = -shiftVector[0];
        shift[1] = -shiftVector[1];
        shift[2] = -shiftVector[2];
        std::shared_ptr<IndicatorTranslate3D<T>> newStlReaderOrigin = std::make_shared<IndicatorTranslate3D<T>>( shift,*newStlReader );

        // Move Domain to Desired Position
        T xPosition = std::stof(_geometryData[1]);
        T yPosition = std::stof(_geometryData[2]);
        T zPosition = std::stof(_geometryData[3]);
      
        shift[0] = xPosition;
        shift[1] = yPosition;
        shift[2] = zPosition;
        std::shared_ptr<IndicatorTranslate3D<T>> newStlReaderShifted = std::make_shared<IndicatorTranslate3D<T>>( shift,*newStlReaderOrigin );
        std::shared_ptr<IndicatorLayer3D<T>> newExtendedDomain = std::make_shared<IndicatorLayer3D<T>>( *newStlReaderOrigin, L );
        
        int noOfCuboids = 1;
        std::shared_ptr<CuboidGeometry3D<T>> newCuboidGeometry = std::make_shared<CuboidGeometry3D<T>>( *newExtendedDomain, L, noOfCuboids, "volume" );
        std::shared_ptr<HeuristicLoadBalancer<T>> newLoadBalancer = std::make_shared<HeuristicLoadBalancer<T>>( *newCuboidGeometry );
        std::shared_ptr<SuperGeometry<T,3>> newSuperGeometry = std::make_shared<SuperGeometry<T,3>>( *newCuboidGeometry, *newLoadBalancer );

        _oldSimulationData[1] = std::to_string(_converterProperties[8]);
        _changedSimulationData[1] = std::to_string(_converterProperties[8]);

        getGeometryData = false;
        while(!getGeometryData)
        {
          getGeometryData = readGeometryData();
          std::this_thread::sleep_for(std::chrono::seconds(2));
        }

        prepareGeometry( *newConverter, *newSuperGeometry, *newExtendedDomain, *newStlReaderShifted, filePaths );

        clout << "Prepare Geometry Successful after Change." << std::endl;

        std::shared_ptr<SuperLattice<T, NSDESCRIPTOR>> NSlatticeNew = std::make_shared<SuperLattice<T, NSDESCRIPTOR>>(*newSuperGeometry);
        std::shared_ptr<SuperLattice<T, TDESCRIPTOR>> ADlatticeNew = std::make_shared<SuperLattice<T, TDESCRIPTOR>>(*newSuperGeometry);

        clout << "Generated Lattices after Change." << std::endl;

        SuperLatticeVelocity3D<T, NSDESCRIPTOR> velocityChange( *NSlattice );
        AnalyticalFfromSuperF3D<T> newfieldVelocity(velocityChange, true);
        SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperatureChange(*ADlattice, *converter);
        AnalyticalFfromSuperF3D<T> newfieldTemperature(temperatureChange, true);

        clout << "Preparing Lattices after Change." << std::endl;
        
        if(std::abs(std::stof(_oldSimulationData[0]) - std::stof(_changedSimulationData[0])) > 0.01)
        {
          // Change in Mesh Quality
          SuperLatticeVelocity3D<T, NSDESCRIPTOR> velocityChangeInitial( *NSlatticeNew );
          AnalyticalFfromSuperF3D<T> newfieldVelocityInitial(velocityChangeInitial, true);
          SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperatureChangeInitial(*ADlatticeNew, *converter);
          AnalyticalFfromSuperF3D<T> newfieldTemperatureInitial(temperatureChangeInitial, true);

          _oldSimulationData[0] = std::to_string(_converterProperties[0]);
          _changedSimulationData[0] = std::to_string(_converterProperties[0]);

          prepareLattice(  *NSlatticeNew, *ADlatticeNew, *newConverter, *newSuperGeometry, newfieldVelocityInitial, newfieldTemperatureInitial);
        }
        else
        {
          clout << "Keeping the old field..." << std::endl;
          prepareLattice(  *NSlatticeNew, *ADlatticeNew, *newConverter, *newSuperGeometry, newfieldVelocity, newfieldTemperature);
        }

        stlReader.swap(newStlReader);
        stlReaderOrigin.swap(newStlReaderOrigin);
        stlReaderShifted.swap(newStlReaderShifted);
        extendedDomain.swap(newExtendedDomain);
        cuboidGeometry.swap(newCuboidGeometry);
        loadBalancer.swap(newLoadBalancer);
        superGeometry.swap(newSuperGeometry);

        NSlattice.swap(NSlatticeNew);
        NSlatticeNew.reset();
        ADlattice.swap(ADlatticeNew);
        ADlatticeNew.reset();

        converter.swap(newConverter);
        newConverter.reset();

        newStlReader.reset();
        newStlReaderOrigin.reset();
        newStlReaderShifted.reset();
        newExtendedDomain.reset();
        newCuboidGeometry.reset();
        newLoadBalancer.reset();
        newSuperGeometry.reset();

        goto label;
      }
    }
  }
  NSlattice->getStatistics().print( converter->getLatticeTime( maxPhysT ),converter->getPhysTime( converter->getLatticeTime( maxPhysT ) ) );

  timer.stop();
  timer.printSummary();
}
