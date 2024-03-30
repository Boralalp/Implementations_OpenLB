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

/* NaturalConvection.cpp:
 * Example Usage of the Heat Source Boundary with the Simple Cabinet Model.
 * This is the simplified version of the Heat Source Boundary used in the Master`s Thesis of Berkay Oralalp.
 * For more accurate results, it is recommended to treat each surface seperatly based on the velocities measured. 
 */

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

using NSDESCRIPTOR = D3Q19<FORCE,TAU_EFF>;
using TDESCRIPTOR = D3Q19<VELOCITY,TAU_EFF>;

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>

// Parameters for the simulation setup
int N = 100;             // resolution of the model
const T cL = .5;         // characteristic length
T L = cL/N;
const T maxPhysT = 500;  // max. simulation time in s, SI unit

const T lengthX = .712;  // Domain x Length [m]
const T lengthY = 2.112; // Domain y Length [m]
const T lengthZ = .512;  // Domain z Length [m]

// Positions of the Objects are initialized as a pointers for 
// run-time objects position change (Not Used for This Example)

auto centerObsX1 = std::make_shared<T>(.6);
auto centerObsY1 = std::make_shared<T>(1.3);
auto centerObsZ1 = std::make_shared<T>(0.);
const T lengthObsX1 = 0.03;
const T lengthObsY1 = 0.08;
const T lengthObsZ1 = 0.2;

auto centerObsX2 = std::make_shared<T>(.5);
auto centerObsY2 = std::make_shared<T>(.4);
auto centerObsZ2 = std::make_shared<T>(0.);
const T lengthObsX2 = 0.04;
const T lengthObsY2 = 0.09;
const T lengthObsZ2 = 0.22;

auto centerObsX3 = std::make_shared<T>(.1);
auto centerObsY3 = std::make_shared<T>(1.3);
auto centerObsZ3 = std::make_shared<T>(0.);
const T lengthObsX3 = 0.05;
const T lengthObsY3 = 0.1;
const T lengthObsZ3 = 0.25;

const T Thot = 317.;        // Max Temperature Estimated
const T Tcold = 300;        // Min Temperature

const T Pr = 0.71;         // Prandtl number

T increment = 0.2;
T smagoConst = 0.1;

// Power of objects in [W]
T power3 = 5;
T power2 = 2.376;
T power1 = 0.96;

// Initial start values for Heat Sources [Kelvin]
auto Tvalue1 = std::make_shared<T>(305);
auto Tvalue2 = std::make_shared<T>(305);
auto Tvalue3 = std::make_shared<T>(305);

// Values to read through the folder
T timerStop;
T prandtlTurbulence = 0.87;
T convergenceTime = 100;
T maxVel = 0.2;

// Stores geometry information in form of material numbers
void prepareGeometry( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                      SuperGeometry<T,3>& superGeometry,
                      std::shared_ptr<IndicatorF3D<T>> Obs1,
                      std::shared_ptr<IndicatorF3D<T>> Obs2,
                      std::shared_ptr<IndicatorF3D<T>> Obs3)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1,1} );
  superGeometry.clean();

  // Set material number for Obstacles
  superGeometry.rename( 1,3, Obs1 );
  superGeometry.rename( 1,4, Obs2 );
  superGeometry.rename( 1,5, Obs3 );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();


  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T,TDESCRIPTOR>& ADlattice,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                     SuperGeometry<T,3>& superGeometry,
                     std::shared_ptr<IndicatorF3D<T>> Obs1,
                     std::shared_ptr<IndicatorF3D<T>> Obs2,
                     std::shared_ptr<IndicatorF3D<T>> Obs3,
                     AnalyticalConst3D<T,T>& uF,
                     AnalyticalConst3D<T,T>& rhoF)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
 
  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();

  AnalyticalConst3D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst3D<T,T> T_hot(converter.getLatticeTemperature(Thot));

  AnalyticalConst3D<T,T> T_hot1(converter.getLatticeTemperature(*Tvalue1));
  AnalyticalConst3D<T,T> T_hot2(converter.getLatticeTemperature(*Tvalue2));
  AnalyticalConst3D<T,T> T_hot3(converter.getLatticeTemperature(*Tvalue3));

  // Material 1 is the Fluid
  // Material 2 is the Wall
  // Material 3 is the Heat Source 1
  // Material 4 is the Heat Source 2
  // Material 5 is the Heat Source 3

  NSlattice.defineDynamics<ExternalTauEffLESForcedBGKdynamics<T,NSDESCRIPTOR>>(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}));
  ADlattice.defineDynamics<ExternalTauEffLESBGKadvectionDiffusionDynamics>(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}));
  NSlattice.setParameter<collision::LES::Smagorinsky>(smagoConst);
  ADlattice.setParameter<collision::LES::Smagorinsky>(smagoConst);
 
  // Define Gravity
  std::vector<T> gravityForce( 3,T() );
  gravityForce[1] = -9.81/converter.getConversionFactorVelocity()*converter.getConversionFactorTime();
  AnalyticalConst3D<T,T> gravitiyForceConst3D( gravityForce );

  // Initialize gravity force
  NSlattice.defineField<FORCE>(superGeometry, 1, gravitiyForceConst3D);
  NSlattice.defineField<FORCE>(superGeometry, 2, gravitiyForceConst3D);

  setLocalVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, NSomega, superGeometry.getMaterialIndicator({2,3,4,5}));

  // ADE lattice definition
  setAdvectionDiffusionTemperatureBoundary<T, TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 2);
  setAdvectionDiffusionTemperatureBoundary<T, TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 3);
  setAdvectionDiffusionTemperatureBoundary<T, TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 4);
  setAdvectionDiffusionTemperatureBoundary<T, TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 5);

  // Initial Values of NSE lattice
  NSlattice.defineRho( superGeometry, 1, rhoF);
  NSlattice.defineU( superGeometry, 1, uF );
  NSlattice.iniEquilibrium( superGeometry, 1, rhoF, uF );

  // Initial Values of ADE lattice
  ADlattice.defineRho(superGeometry, 1, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 1, T_cold, uF);
  ADlattice.defineRho(superGeometry, 2, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 2, T_cold, uF);
  ADlattice.defineRho(superGeometry, 3, T_hot1);
  ADlattice.iniEquilibrium(superGeometry, 3, T_hot1, uF);
  ADlattice.defineRho(superGeometry, 4, T_hot2);
  ADlattice.iniEquilibrium(superGeometry, 4, T_hot2, uF);
  ADlattice.defineRho(superGeometry, 5, T_hot3);
  ADlattice.iniEquilibrium(superGeometry, 5, T_hot3, uF);

  //TAU_EFF Smagorinsky Parameters
  AnalyticalConst3D<T,T> tauNS(1./NSomega);
  AnalyticalConst3D<T,T> tauAD(1./Tomega);

  NSlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}), tauNS );
  ADlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}), tauAD );

  NSlattice.setParameter<descriptors::OMEGA>(NSomega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);

  /// Make the lattice ready for simulation
  clout << "initialize ..." << std::endl;
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}


// Compute Heat Source
void computePowerInOne(SuperGeometry<T,3>& superGeometry,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice,
                 int materialNumber,
                 std::shared_ptr<IndicatorF3D<T>> Object,
                 T power,
                 std::shared_ptr<T> Tvalue,
                 ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter)
{
  ADlattice.setProcessingContext(ProcessingContext::Evaluation);
  NSlattice.setProcessingContext(ProcessingContext::Evaluation);

  int voxel = 0, material = 0;
  T T_x = 0, T_x1 = 0, T_x2 = 0;
  T q = 0;

  std::vector<T> velocityAroundCells = {0,0,0};

  T Nx = 0,Ny = 0, Nz = 0;

  Vector<T,3> startCorner;
  Vector<T,3> endCorner;
  startCorner = Object->getMin();
  endCorner = Object->getMax();

  Nx = (endCorner[0] - startCorner[0]) / L;
  Ny = (endCorner[1] - startCorner[1]) / L;
  Nz = (endCorner[2] - startCorner[2]) / L;

  std::cout << "The Dimensions are (" << Nx << ", " << Ny << ", " << Nz << ")" << std::endl;

  //Loop Around the Object
  //The loop can be optimized by changing the increment to the length of the object
  for (int iC = 0; iC < NSlattice.getLoadBalancer().size(); iC++) 
  {
    for (int iZ = startCorner[2]/L - 5; iZ < endCorner[2]/L + 5; ++iZ)
    {
      for (int iY = startCorner[1]/L - 5; iY < endCorner[1]/L + 5; ++iY) 
      {
        for (int iX = startCorner[0]/L - 5; iX <= endCorner[0]/L + 5; ++iX)
        {
          material = superGeometry.getBlockGeometry(iC).getMaterial(iX,iY,iZ);

          if ( material == materialNumber ) 
          {
            // Discrete Normals can be used. However, it has been seen that it is hard for the mesher 
            // to obtain normals to boundaries for complex geometries.
            // For Natural Convection, assumptions can be made on the direction of the flow.
            if(superGeometry.getBlockGeometry(iC).getMaterial(iX - 1, iY, iZ) == 1)
            {
              /*

              * Example code on getting cell velocities to get more accurate results on the Nusselts number.
              
              const auto extVel = ADlattice.getBlock(iC).get(iX - 1, iY, iZ).template getField<descriptors::VELOCITY>();
              velocityAroundCells[0] = std::abs(extVel[0]);
              velocityAroundCells[1] = std::abs(extVel[1]);
              velocityAroundCells[2] = std::abs(extVel[2]);

              */

              T_x = ADlattice.getBlock(iC).get(iX - 1,iY,iZ).computeRho();
              T_x1 = ADlattice.getBlock(iC).get(iX - 2, iY, iZ).computeRho();
              T_x2 = ADlattice.getBlock(iC).get(iX - 3, iY, iZ).computeRho();
              q += ((3.0*T_x - 4.0*T_x1 + 1.0*T_x2)/2.0)*Ny;
              voxel++;
            }
            else if(superGeometry.getBlockGeometry(iC).getMaterial(iX + 1, iY, iZ) == 1)
            {
              T_x = ADlattice.getBlock(iC).get(iX + 1,iY,iZ).computeRho();
              T_x1 = ADlattice.getBlock(iC).get(iX + 2, iY, iZ).computeRho();
              T_x2 = ADlattice.getBlock(iC).get(iX + 3, iY, iZ).computeRho();
              
              q += (3.0*T_x - 4.0*T_x1 + 1.0*T_x2)/2.0*Ny;
              voxel++;
            }

            else if(superGeometry.getBlockGeometry(iC).getMaterial(iX, iY - 1, iZ) == 1)
            {
              T_x = ADlattice.getBlock(iC).get(iX, iY - 1, iZ).computeRho();
              T_x1 = ADlattice.getBlock(iC).get(iX, iY - 2, iZ).computeRho();
              T_x2 = ADlattice.getBlock(iC).get(iX, iY - 3, iZ).computeRho();
              
              q += (3.0*T_x - 4.0*T_x1 + 1.0*T_x2)/2.0*Nz;
              voxel++;
            }
            else if(superGeometry.getBlockGeometry(iC).getMaterial(iX, iY + 1, iZ) == 1)
            {
              T_x = ADlattice.getBlock(iC).get(iX, iY + 1, iZ).computeRho();
              T_x1 = ADlattice.getBlock(iC).get(iX, iY + 2, iZ).computeRho();
              T_x2 = ADlattice.getBlock(iC).get(iX, iY + 3, iZ).computeRho();
              
              q += (3.0*T_x - 4.0*T_x1 + 1.0*T_x2)/2.0*Nz;
              voxel++;
            }
            else if(superGeometry.getBlockGeometry(iC).getMaterial(iX, iY, iZ - 1) == 1)
            {
              T_x = ADlattice.getBlock(iC).get(iX, iY, iZ - 1).computeRho();
              T_x1 = ADlattice.getBlock(iC).get(iX, iY, iZ - 2).computeRho();
              T_x2 = ADlattice.getBlock(iC).get(iX, iY, iZ - 3).computeRho();
              
              q += (3.0*T_x - 4.0*T_x1 + 1.0*T_x2)/2.0*Ny;
              voxel++;
            }
            else if(superGeometry.getBlockGeometry(iC).getMaterial(iX, iY, iZ + 1) == 1)
            {
              T_x = ADlattice.getBlock(iC).get(iX, iY, iZ + 1).computeRho();
              T_x1 = ADlattice.getBlock(iC).get(iX, iY, iZ + 2).computeRho();
              T_x2 = ADlattice.getBlock(iC).get(iX, iY, iZ + 3).computeRho();
              
              q += (3.0*T_x - 4.0*T_x1 + 1.0*T_x2)/2.0*Ny;
              voxel++;
            }
            else
            {
              // Inside the object! Nothing to do...
            }
          }
        }
      }
    }
  }

  T nusseltNumber = q / (T)voxel;

  std::cout << "Nusselt number value for material " << materialNumber << " is " << nusseltNumber << std::endl;

  T faceArea = (2 * (Nx * Ny) + 2 * (Ny * Nz) + 2 * (Nx * Nz) ) * std::pow(L, 2.);
  T volume = (Nx * Ny * Nz) * std::pow(L, 3.);
  T charLength = std::pow(volume, 1/3.);

  std::cout << "Face area for material " << materialNumber << " is " << faceArea << " m^2." << std::endl;
  std::cout << "charLength for material " << materialNumber << " is " << charLength << " m." << std::endl;
  T currentPowerOutput = (*Tvalue - converter.getCharPhysLowTemperature()) * (0.025684/charLength) * faceArea * nusseltNumber;

  std::cout << "Goal PowerOutput for material " << materialNumber << " is " << power << " Watts." << std::endl;
  std::cout << "currentPowerOutput for material " << materialNumber << " is " << currentPowerOutput << " Watts." << std::endl;
  std::cout << "Previous Temperature value for material " << materialNumber << " is " << *Tvalue << " Kelvin." << std::endl;

  if(power - currentPowerOutput < 0.5)
  {
    increment = 0.05;
  }
  else
  {
    increment = 1;
  }

  if(power > currentPowerOutput)
  {
    *Tvalue += increment;
  }
  else
  {
    *Tvalue -= increment;
  }

  std::cout << "New Temperature value for material " << materialNumber << " is " << *Tvalue << " Kelvin." << std::endl;
  std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

  // Update ADE Lattice
  AnalyticalConst3D<T,T> T_value(converter.getLatticeTemperature(*Tvalue));
  ADlattice.defineRho(superGeometry, materialNumber, T_value);

  ADlattice.setProcessingContext(ProcessingContext::Simulation);
  NSlattice.setProcessingContext(ProcessingContext::Simulation);
}


// Computes the pressure drop between the voxels before and after the cylinder
void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 std::shared_ptr<IndicatorF3D<T>> Obs1,
                 std::shared_ptr<IndicatorF3D<T>> Obs2,
                 std::shared_ptr<IndicatorF3D<T>> Obs3)
{
  OstreamManager clout( std::cout,"getResults" );
 
  SuperVTMwriter3D<T> vtmWriter( "simpleCabinet" );
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

  if ( iT==0 ) 
  {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry( NSlattice, superGeometry );
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid( NSlattice );
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank( NSlattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
 
    vtmWriter.createMasterFile();
  }

  // Computes new temperature values and prints the results
  if ( iT%vtkIter == 0 && iT >= converter.getLatticeTime(0.)) 
  {
    {
      computePowerInOne(superGeometry, NSlattice, ADlattice, 3, Obs1, power1, Tvalue1, converter);
      computePowerInOne(superGeometry, NSlattice, ADlattice, 4, Obs2, power2, Tvalue2, converter);
      computePowerInOne(superGeometry, NSlattice, ADlattice, 5, Obs3, power3, Tvalue3, converter);

      SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( NSlattice, converter );
      SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);

      vtmWriter.addFunctor( velocity );
      vtmWriter.addFunctor( temperature );
      vtmWriter.write( iT );
    }
  }
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  //Read Some Parameters from a .txt file.
  std::ifstream inputFile("input.txt");

  if (inputFile.is_open()) 
  {
      std::string line;
      if (std::getline(inputFile, line))
      {
          std::istringstream iss(line);

          std::string intString;
          if (std::getline(iss, intString, ',')) {
              N = std::stoi(intString);
          } else {
              std::cerr << "Error reading integer value." << std::endl;
              return 1;
          }

          std::string floatString;
          if (std::getline(iss, floatString, ',')) {
              smagoConst = std::stof(floatString);
          } else {
              std::cerr << "Error reading float value." << std::endl;
              return 1;
          }

          std::string floatString2;
          if (std::getline(iss, floatString2, ',')) {
              prandtlTurbulence = std::stof(floatString2);
          } else {
              std::cerr << "Error reading float value." << std::endl;
              return 1;
          }

          std::string floatString3;
          if (std::getline(iss, floatString3, ',')) {
              convergenceTime = std::stof(floatString3);
          } else {
              std::cerr << "Error reading float value." << std::endl;
              return 1;
          }

          std::string floatString4;
          if (std::getline(iss, floatString4)) {
              maxVel = std::stof(floatString4);
          } else {
              std::cerr << "Error reading float value." << std::endl;
              return 1;
          }
      }
      else 
      {
          std::cerr << "Error reading line from the file." << std::endl;
      }

      // Close the file
      inputFile.close();
  }

  //Reading opeartion done.

  L = cL/N;

  //Converter for ADE and NSE lattices
  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
    (T) cL/N,                   // physDeltaX
    (T) 2.*0.05/(maxVel)*cL/N,  // physDeltaT
    (T) cL,                     // charPhysLength
    (T) maxVel,                 // charPhysVelocity
    (T) 15.126e-6,              // Viscosity
    (T) 1.0,                    // Density
    (T) 25.684e-3,
    (T) Pr * 25.684e-3 / 15.126e-6 / 1.0,
    (T) 0.00341,
    (T) Tcold,                  // Min Temperature
    (T) Thot                    // Max Temperature
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("NaturalConvection3d");
 
  // Prepare Domain Geometry
  Vector<T,3> extend( lengthX,lengthY,lengthZ );
  Vector<T,3> origin;
  IndicatorCuboid3D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  clout << "Number Of Cuboids Used: " << noOfCuboids << std::endl;
  CuboidGeometry3D<T> cuboidGeometry( cuboid, L, noOfCuboids, "volume" );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  // Prepare Heat Source Geometry
  Vector<T,3> originObj1( *centerObsX1,*centerObsY1,*centerObsZ1);
  std::shared_ptr<Vector<T,3>> extendObj1 = std::make_shared<Vector<T,3>>( lengthObsX1,lengthObsY1,lengthObsZ1 );
  std::shared_ptr<IndicatorF3D<T>> Obs1 = std::make_shared<IndicatorCuboid3D<T>>( *extendObj1, originObj1 );

  Vector<T,3> originObj2( *centerObsX2,*centerObsY2,*centerObsZ2);
  std::shared_ptr<Vector<T,3>> extendObj2 = std::make_shared<Vector<T,3>>( lengthObsX2,lengthObsY2,lengthObsZ2 );
  std::shared_ptr<IndicatorF3D<T>> Obs2 = std::make_shared<IndicatorCuboid3D<T>>( *extendObj2, originObj2 );

  Vector<T,3> originObj3( *centerObsX3,*centerObsY3,*centerObsZ3);
  std::shared_ptr<Vector<T,3>> extendObj3 = std::make_shared<Vector<T,3>>( lengthObsX3,lengthObsY3,lengthObsZ3 );
  std::shared_ptr<IndicatorF3D<T>> Obs3 = std::make_shared<IndicatorCuboid3D<T>>( *extendObj3, originObj3 );

  // Initial conditions
  auto rhoF = std::make_shared<AnalyticalConst3D<T,T>>( 1 );
  std::vector<T> velocityZero( 3,T( 0 ) );
  auto uF = std::make_shared<AnalyticalConst3D<T,T>>( velocityZero );
  
  // Prepare Geometry
  prepareGeometry( converter, superGeometry, Obs1, Obs2, Obs3 );

  // Prepare Lattice
  SuperLattice<T, NSDESCRIPTOR> NSlattice( superGeometry );
  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);

  prepareLattice(  NSlattice, ADlattice, converter, superGeometry, Obs1, Obs2, Obs3, *uF, *rhoF);

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() * converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

  const T preFactor = smagoConst*smagoConst
                    * descriptors::invCs2<T,NSDESCRIPTOR>()*descriptors::invCs2<T,NSDESCRIPTOR>()
                    * 2*util::sqrt(2);
                
  // Couple Lattices
  SuperLatticeCoupling coupling(
    SmagorinskyBoussinesqCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{},  ADlattice);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::T0>(
    converter.getLatticeTemperature(Tcold));
  coupling.setParameter<SmagorinskyBoussinesqCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,3>{0.0,1.0,0.0});
  coupling.setParameter<SmagorinskyBoussinesqCoupling::SMAGORINSKY_PREFACTOR>(preFactor);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::PR_TURB>(prandtlTurbulence);
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_NSE>(
    converter.getLatticeRelaxationFrequency());
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_ADE>(
    converter.getLatticeThermalRelaxationFrequency());

  // Main Loop with Timer
  clout << "Starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) 
  {
    // Collide and Stream Execution
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();

    coupling.execute();

    // Computation and Output of the Results
    getResults( converter, NSlattice, ADlattice, iT, superGeometry, timer, Obs1, Obs2, Obs3 );
  }

  NSlattice.getStatistics().print( converter.getLatticeTime( maxPhysT ),converter.getPhysTime( converter.getLatticeTime( maxPhysT ) ) );

  timer.stop();
  timer.printSummary();
  
}
