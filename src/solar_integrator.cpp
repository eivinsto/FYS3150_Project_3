#include "integrator.hpp"

#include <string>
#include <armadillo>

integrator::integrator(double dt, std::string integrator) {
  m_dt = dt;
  m_integrator = integrator;
  m_dt_2 = dt/2;     // For use in velocity verlet
  m_dt_sqr_2 = dt*dt/2;  // For use in velocity verlet
}

void integrator::integrateOneStep(SolarSystem& system) {
  if (m_integrator=="Euler") {
    Euler(vec& x, vec& v, vec& a);
  }

  if (m_integrator=="VelocityVerlet") {
    VelocityVerlet(vec& x, vec& v, vec& a);
  }
}

void integrator::Euler(SolarSystem& system) {
  system.calculateForcesAndEnergy();

  for (CelestialBody &body: system.bodies()) {
    body.position += body.velocity*m_dt;
    body.velocity += m_dt * body.force / body.mass;
  }
}

void integrator::VelocityVerlet(SolarSystem& system) {
  system.calculateForcesAndEnergy();
  arma::mat prev_acceleration(system.numberOfBodies(),3);
  int i = 0;

  for (CelestialBody &body: system.bodies()) {
    prev_acceleration.row(i) &= body.force / body.mass;
    body.position += m_dt * body.velocity + m_dt_sqr_2 * prev_acceleration;
    ++i;
  }

  system.calculateForcesAndEnergy()
  i = 0;
  for (CelestialBody &body: system.bodies()) {
    body.velocity += m_dt_2*( prev_acceleration.row(i) + body.force / body.mass)
    ++i;
  }
}
