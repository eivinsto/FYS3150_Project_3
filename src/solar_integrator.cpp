#include "solar_integrator.hpp"

#include <string>
#include <armadillo>

solar_integrator::solar_integrator(double dt, std::string integrator) {
  m_dt = dt;
  m_integrator = integrator;
  m_dt_2 = dt/2;     // For use in velocity verlet
  m_dt_sqr_2 = dt*dt/2;  // For use in velocity verlet
}

void solar_integrator::integrateOneStep(SolarSystem& system) {
  if (m_integrator=="Euler") {
    Euler(system);
  }

  if (m_integrator=="VelocityVerlet") {
    VelocityVerlet(system);
  }
}

void solar_integrator::Euler(SolarSystem& system) {
  system.calculateForcesAndEnergy();

  for (CelestialBody &body: system.bodies()) {
    body.position += body.velocity*m_dt;
    body.velocity += m_dt * body.force / body.mass;
  }
}

void solar_integrator::VelocityVerlet(SolarSystem& system) {
  system.calculateForcesAndEnergy();
  arma::mat prev_acceleration(system.numberOfBodies(),3);
  int i = 0;

  for (CelestialBody &body: system.bodies()) {
    prev_acceleration.row(i) = body.force / body.mass;
    body.position += m_dt * body.velocity + m_dt_sqr_2 * prev_acceleration;
    ++i;
  }

  system.calculateForcesAndEnergy();
  i = 0;
  for (CelestialBody &body: system.bodies()) {
    body.velocity += m_dt_2*( prev_acceleration.row(i) + body.force / body.mass);
    ++i;
  }
}
