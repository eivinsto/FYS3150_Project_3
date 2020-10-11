#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP

#include "celestial_body.hpp"
#include <vector>
#include <armadillo>

class SolarSystem {
public:
  SolarSystem();
  CelestialBody& createCelestialBody(arma::vec& pos, arma::vec& vel, double m);
  int numberOfBodies() const;
  void calculateForcesAndEnergy();
  std::vector<CelestialBody> bodies();

private:
  std::vector<CelestialBody> m_bodies;
  double m_potential_energy;
  double m_kinetic_energy;
};

#endif
