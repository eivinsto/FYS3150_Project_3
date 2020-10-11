#ifndef SOLAR_INTEGRATOR.HPP
#define SOLAR_INTEGRATOR.HPP

#include <string>
#include <armadillo>

class integrator {
public:
  double m_dt;
  double m_dt_sqr_2;
  double m_dt_2;
  std::string m_integrator;
  void integrateOneStep(SolarSystem& system);
  integrator(double dt, std::string integrator);

private:
  void Euler(SolarSystem& system);
  void VelocityVerlet(SolarSystem& system);
};

#endif
