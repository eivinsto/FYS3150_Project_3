#define _USE_MATH_DEFINES
#include "solar_integrator.hpp"
#include "solar_system.hpp"
#include "celestial_body.hpp"
#include <string>
#include <armadillo>
#include <cmath>

using namespace std;

int main(int numArguments, char **arguments) {
  std::string runflag = arguments[1];
  if (runflag == "se") {
    int numTimesteps = 1000;
    if(numArguments >= 2) numTimesteps = atoi(arguments[2]);

    // Reading initial state from file
    SolarSystem solarSystem("../data/earth-sun-init.txt");
    solarSystem.moveToCOFMFrame();

    double dt = 0.001;
    if(numArguments >= 3) dt = atof(arguments[3]);

    solar_integrator integrator(dt, "VelocityVerlet");
    solarSystem.initiateDataFile("../data/positions.xyz", "../data/energies.dat");
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      integrator.integrateOneStep(solarSystem);
      solarSystem.writeToFile();
      solarSystem.writeEnergyToFile();
    }
  }
  return 0;
}
