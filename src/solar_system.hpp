#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP

#include "celestial_body.hpp"
#include <vector>
#include <armadillo>

class SolarSystem {
public:
  SolarSystem();
  SolarSystem(double beta);
  SolarSystem(std::string input_filename);
  SolarSystem(std::string input_filename, double beta);
  CelestialBody& createCelestialBody(arma::vec& pos, arma::vec& vel, double m);
  int numberOfBodies() const;
  void calculateForcesAndEnergy();
  void calculateForcesAndEnergyWithRelativisticCorrection();
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
  double m_beta;
  const double m_c = 63239.7263;  // AU/y
  const double m_rel_constant = 3/(m_c*m_c); // For use in calculateForcesAndEnergyWithRelativisticCorrection
};

#endif
