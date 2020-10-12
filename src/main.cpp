#define _USE_MATH_DEFINES
#include "solar_integrator.hpp"
#include "solar_system.hpp"
#include "celestial_body.hpp"
#include <armadillo>
#include <cmath>

using namespace std;

int main(int numArguments, char **arguments) {
  int numTimesteps = 1000;
  if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

  SolarSystem solarSystem;
  // We create new bodies like this. Note that the createCelestialBody function returns a reference to the newly created body
  // This can then be used to modify properties or print properties of the body if desired
  // Use with: solarSystem.createCelestialBody( position, velocity, mass );

  arma::vec suninit = {0,0,0};
  CelestialBody &sun = solarSystem.createCelestialBody( suninit, suninit, 1.0 );

  // We don't need to store the reference, but just call the function without a left hand side
  arma::vec earthposinit = {1,0,0};
  arma::vec earthvelinit = {0, 2*M_PI, 0};
  solarSystem.createCelestialBody( earthposinit, earthvelinit, 3e-6 );


  double dt = 0.001;
  solar_integrator integrator(dt, "VelocityVerlet");
  solarSystem.initiateDataFile("../data/positions.xyz");
  for(int timestep=0; timestep<numTimesteps; timestep++) {
    integrator.integrateOneStep(solarSystem);
    solarSystem.writeToFile();
  }
  return 0;
}
