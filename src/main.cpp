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
  int write_limit = 1;                              // How often data should be written to file
  string integration_method = "VelocityVerlet";     // Which integration method to use
  string init_file = "../data/sun-earth.init";      // File containing initial conditions
  string positions_file = "../data/positions.xyz";  // Which file to store output positions in
  string energies_file = "../data/energies.dat";    // Which file to store energies and angular momentum output in
  string correction = "nonrel";                     // whether to use relativistic correction of gravity
  double beta = 2.0;                                // exponent for r in gravitational force
  string runflag = "se";                            // runflag for whether to run perihelion precession simulation

  // Reading parameters from command line
  if(numArguments >= 2) numTimesteps = atoi(arguments[1]);
  if(numArguments >= 3) dt = atof(arguments[2]);
  if(numArguments >= 4) write_limit = atoi(arguments[3]);
  if(numArguments >= 5) integration_method = arguments[4];
  if(numArguments >= 6) init_file = arguments[5];
  if(numArguments >= 7) positions_file = arguments[6];
  if(numArguments >= 8) energies_file = arguments[7];
  if(numArguments >= 9) correction = arguments[8];
  if(numArguments >= 10) beta = atof(arguments[9]);
  if(numArguments >= 11) runflag = arguments[10];

  // printing exponent beta
  std::cout << "Current beta = " << beta << std::endl;

  // Reading initial state from file
  SolarSystem solarSystem(init_file, beta);

  // Initiate integrator
  solar_integrator integrator(dt, integration_method);
  solarSystem.initiateDataFile(positions_file, energies_file);

  // moves system to center of mass frame of reference unless relativistic correction is used:
  if (correction == "nonrel") solarSystem.moveToCOFMFrame();
  // precalculates angular momentum for perihelion precession simulation with relativistic correction:
  if (correction == "rel") solarSystem.calAngMom();

  if (runflag == "sm") {
    // runs simulation of sun-mercury system and prints final perihelion to position file,
    double rmin;
    double r;
    arma::vec perihelion;
    CelestialBody *mercury = &(solarSystem.bodies()[1]);

    solarSystem.writeToFile();
    rmin = arma::norm(mercury->position);
    r = rmin;

    // integrates orbit for first 99 years of simulation (assuming numTimesteps*dt = 100 years):
    for(int timestep=0; timestep<numTimesteps - numTimesteps/100; timestep++) {
      integrator.integrateOneStep(solarSystem, correction);
      r = arma::norm(mercury->position);
    }
    rmin = r; // rmin for finding perihelion.

    // integrates remaining year of simulation while looking for perihelion:
    for(int timestep=0; timestep<numTimesteps/100; timestep++) {
      integrator.integrateOneStep(solarSystem, correction);
      r = arma::norm(mercury->position);

      if (r < rmin) {
        // storing smallest distance between Mercury and Sun:
        perihelion = mercury->position;
        rmin = r;
      }
    }
    // writing perihelion to position data file:
    mercury->position = perihelion;
    solarSystem.writeToFile();

  } else {
    // running all other simulations:
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      integrator.integrateOneStep(solarSystem, correction);

      // allowing for limited writes to file to speed up integration:
      if (timestep%write_limit == 0){
        solarSystem.writeToFile();
        solarSystem.writeEnergyToFile();
      }
    }
  }

  return 0;
}
