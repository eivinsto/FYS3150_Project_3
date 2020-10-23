#include "solar_integrator.hpp"

#include <string>
#include <armadillo>

/**
* Constructor that takes a timestep dt and a string that tells the class
* which integrator to use. Currently the Euler method and the velocity
* verlet method are implemented.
*
* @dt Timestep.
* @integrator String that tells the class which integrator to use. Can be either
* "Euler" or "VelocityVerlet".
*/
solar_integrator::solar_integrator(double dt, std::string integrator) {
  m_dt = dt;
  m_integrator = integrator;
  m_dt_2 = dt/2;              // For use in velocity verlet
  m_dt_sqr_2 = dt*dt/2;       // For use in velocity verlet
}

/**
* Function that integrates one step in the SolarSystem object specified.
* Chooses which method used depending on keyword specified in constructor.
*
* @system SolarSystem object to integrate.
* @relOrNonRel string containing whether to use relativistic corrections.
*/
void solar_integrator::integrateOneStep(SolarSystem& system, std::string relOrNonRel) {
  // If-test determining which integration method to use
  if (m_integrator=="Euler") {
    Euler(system, relOrNonRel);
  }

  if (m_integrator=="VelocityVerlet") {
    VelocityVerlet(system, relOrNonRel);
  }
}

/**
* Function that calls integrateOneStep without relativistic corrections.
* Chooses which method used depending on keyword specified in constructor.
*
* @system SolarSystem object to integrate.
*/
void solar_integrator::integrateOneStep(SolarSystem& system) {
  solar_integrator::integrateOneStep(system, "nonrel");
}

/**
* Function that integrates one step with the Euler method. This function is
* by the integrateOneStep() function internally if the method to be used was
* specified to be the Euler method in the constructor.
*
* @system SolarSystem object to be integrated.
* @relOrNonRel string containing flag for whether to use relativistic correction.
*/
void solar_integrator::Euler(SolarSystem& system, std::string relOrNonRel) {
  // Calculate forces and energy
  if (relOrNonRel == "rel") {
    system.calculateForcesWithRelativisticCorrection();
  } else if (relOrNonRel == "nonrel") {
    system.calculateForcesAndEnergy();
  }

  // Integrate every body one step
  for (CelestialBody &body: system.bodies()) {
    body.position += body.velocity*m_dt;
    body.velocity += m_dt * body.force / body.mass;
  }
}

// /**
// * Function calling Euler with no relativistic correction.
// *
// * @system SolarSystem object to be integrated.
// */
// void solar_integrator::Euler(SolarSystem& system) {
//   solar_integrator::Euler(system, "nonrel");
// }

/**
* Function that integrates one step with the velocity verlet method. This
* function is by the integrateOneStep() function internally if the method to be
* used was specified to be the velocity verlet method in the constructor.
*
* @system SolarSystem object to be integrated.
*/
void solar_integrator::VelocityVerlet(SolarSystem& system, std::string relOrNonRel) {
  // Calculate forces and energy
  CelestialBody body1 = system.bodies()[1];
  if (relOrNonRel == "rel") {
    system.calculateForcesWithRelativisticCorrection();
  } else if (relOrNonRel == "nonrel") {
    system.calculateForcesAndEnergy();
  }

  // Initializing matrix to store the acceleration in current timestep
  arma::mat prev_acceleration(3,system.numberOfBodies());

  // Integrates the position of the bodies one step
  int i = 0;
  for (CelestialBody &body: system.bodies()) {
    prev_acceleration.col(i) = body.force / body.mass;
    body.position += m_dt * body.velocity + m_dt_sqr_2 * prev_acceleration.col(i);
    ++i;
  }

  // Calculate forces and energy again (to get acceleration in next timestep)
  if (relOrNonRel == "rel") {
    system.calculateForcesWithRelativisticCorrection();
  } else if (relOrNonRel == "nonrel") {
    system.calculateForcesAndEnergy();
  }

  // Integrates the velocity of the bodies one step using the acceleration in the current and the next timestep
  i = 0;
  for (CelestialBody &body: system.bodies()) {
    body.velocity += m_dt_2*( prev_acceleration.col(i) + body.force / body.mass);
    ++i;
  }
}

// /**
// * Function calling VelocityVerlet with no relativistic correction.
// *
// * @system SolarSystem object to be integrated.
// */
// void solar_integrator::VelocityVerlet(SolarSystem& system) {
//   solar_integrator::VelocityVerlet(system, "nonrel");
// }
