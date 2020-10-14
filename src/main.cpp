#define _USE_MATH_DEFINES
#include "solar_integrator.hpp"
#include "solar_system.hpp"
#include "celestial_body.hpp"
#include <string>
#include <armadillo>
#include <cmath>

using namespace std;

int main(int numArguments, char **arguments) {
  int numTimesteps = 1000;                          // Amount of timesteps
  double dt = 0.001;                                // Timestep
  string init_file = "../data/earth-sun-init.txt";  // File containing initial conditions
  int write_limit = 1;                              // How often data should be written to file
  string integration_method = "VelocityVerlet";     // Which integration method to use
  string positions_file = "../data/positions.xyz";  // Which file to store output positions in
  string energies_file = "../data/energies.dat";    // Which file to store energies and angular momentum output in

  // Reading variables from command line
  if(numArguments >= 2) numTimesteps = atoi(arguments[1]);
  if(numArguments >= 3) dt = atod(arguments[2]);
  if(numArguments >= 4) init_file = arguments[3];
  if(numArguments >= 5) write_limit = atoi(arguments[4]);
  if(numArguments >= 6) integration_method = arguments[5];
  if(numArguments >= 7) positions_file = "../data/" + arguments[6];
  if(numArguments >= 8) energies_file = "../data/" + arguments[7];

  // Reading initial state from file
  SolarSystem solarSystem(init_file);

  // Initiate integrator
  solar_integrator integrator(dt, integration_method);
  solarSystem.initiateDataFile("../data/positions.xyz", "../data/energies.dat");
  for(int timestep=0; timestep<numTimesteps; timestep++) {
    integrator.integrateOneStep(solarSystem);

    if (i%write_limit == 0){
    solarSystem.writeToFile();
    solarSystem.writeEnergyToFile();
    }
  }
  return 0;
}
