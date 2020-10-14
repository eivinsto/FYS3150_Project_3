#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP

#define _USE_MATH_DEFINES
#include "celestial_body.hpp"
#include <cmath>
#include <vector>
#include <armadillo>

class SolarSystem {
public:
  // Public methods
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
  // Private members
  std::vector<CelestialBody> m_bodies;       // Vector storing the celestial body objects
  double m_potential_energy;                 // Total potential energy ([mass of the sun] AU^2 y^-2)
  double m_kinetic_energy;                   // Total kinetic energy ([mass of the sun] AU^2 y^-2)
  arma::vec m_angular_momentum;              // Total angular momentum vector ([mass of the sun] AU^2 y^-1)
  std::string m_positions_filename;          // Name of file to write positions to
  std::ofstream m_file;                      // File to write positions to
  std::string m_energy_filename;             // Name of file to write energy and angular momentum to
  std::ofstream m_energy_file;               // File to write energy and angular momentum to
  double m_beta;                             // Exponent in force
  const double m_G = 4*M_PI*M_PI;            // Gravitational constant (AU^3 y^-2 [mass of the sun]^-1)
  const double m_c = 63239.7263;             // Speed of light (AU y^-1)
  const double m_rel_constant = 3/(m_c*m_c); // For use in calculateForcesAndEnergyWithRelativisticCorrection (y^2 AU^-2)
};

#endif
