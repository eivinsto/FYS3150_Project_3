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
  string init_file = "../data/sun-earth.init";  // File containing initial conditions
  string positions_file = "../data/positions.xyz";  // Which file to store output positions in
  string energies_file = "../data/energies.dat";    // Which file to store energies and angular momentum output in
  string correction = "nonrel";
  double beta = 2.0;
  string runflag = "se";

  // Reading variables from command line
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
  std::cout << "Current beta = " << beta << std::endl;
  // Reading initial state from file
  SolarSystem solarSystem(init_file, beta);

  // Initiate integrator
  solar_integrator integrator(dt, integration_method);
  solarSystem.initiateDataFile(positions_file, energies_file);

  if (correction == "nonrel") solarSystem.moveToCOFMFrame();
  if (correction == "rel") solarSystem.calAngMom();

  if (runflag == "sm") {
    double rmin;
    double r;
    arma::vec perihelion;
    CelestialBody *mercury = &(solarSystem.bodies()[1]);

    solarSystem.writeToFile();
    rmin = arma::norm(mercury->position);


    for(int timestep=0; timestep<numTimesteps - numTimesteps/100; timestep++) {
      integrator.integrateOneStep(solarSystem, correction);
      r = arma::norm(mercury->position);
    }
    rmin = r;

    for(int timestep=0; timestep<numTimesteps/100; timestep++) {
      integrator.integrateOneStep(solarSystem, correction);
      r = arma::norm(mercury->position);

      if (r < rmin) {
        perihelion = mercury->position;
      }
    }
    mercury->position = perihelion;
    solarSystem.writeToFile();

  } else {
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      integrator.integrateOneStep(solarSystem, correction);

      if (timestep%write_limit == 0){
        solarSystem.writeToFile();
        solarSystem.writeEnergyToFile();
      }
    }
  }
  return 0;
}
