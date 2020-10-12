#include "solar_integrator.hpp"
#include "solar_system.hpp"
#include "celestial_body.hpp"
#include <armadillo>

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

    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(int i = 0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
    }

    double dt = 0.001;
    solar_integrator integrator(dt, "VelocityVerlet");
    for(int timestep=0; timestep<numTimesteps; timestep++) {
        integrator.integrateOneStep(solarSystem);
        solarSystem.writeToFile("../data/positions.xyz");
    }

    cout << "I just created my first solar system that has " << solarSystem.bodies().size() << " objects." << endl;
    return 0;
  return 0;
}
