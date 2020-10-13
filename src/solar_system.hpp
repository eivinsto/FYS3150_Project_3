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
  std::vector<CelestialBody> &bodies();
  double potentialEnergy();
  double kineticEnergy();
  arma::vec angularMomentum();
  void initiateDataFile(std::string positions_filename);
  void initiateDataFile(std::string positions_filename, std::string energy_filename);
  void writeToFile();
  void writeEnergyToFile();

private:
  std::vector<CelestialBody> m_bodies;
  double m_potential_energy;
  double m_kinetic_energy;
  arma::vec m_angular_momentum;
  std::string m_positions_filename;
  std::ofstream m_file;
  std::string m_energy_filename;
  std::ofstream m_energy_file;
};

#endif
