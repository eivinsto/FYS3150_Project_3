#ifndef SOLAR_INTEGRATOR_HPP
#define SOLAR_INTEGRATOR_HPP

#include <string>
#include <armadillo>
#include "solar_system.hpp"

class solar_integrator {
public:
  double m_dt;
  double m_dt_sqr_2;
  double m_dt_2;
  std::string m_integrator;
  void integrateOneStep(SolarSystem& system);
  solar_integrator(double dt, std::string integrator);

private:
  void Euler(SolarSystem& system);
  void VelocityVerlet(SolarSystem& system);
};

#endif
