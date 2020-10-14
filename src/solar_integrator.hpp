#ifndef SOLAR_INTEGRATOR_HPP
#define SOLAR_INTEGRATOR_HPP

#include <string>
#include <armadillo>
#include "solar_system.hpp"

class solar_integrator {
public:
  // Public function
  void integrateOneStep(SolarSystem& system);

  // Constructor
  solar_integrator(double dt, std::string integrator);

private:
  // Private variables
  double m_dt;
  std::string m_integrator;
  double m_dt_sqr_2;
  double m_dt_2;

  // Private functions
  void Euler(SolarSystem& system);
  void VelocityVerlet(SolarSystem& system);
};

#endif
