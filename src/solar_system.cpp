#define _USE_MATH_DEFINES

#include <cmath>
#include "solar_system.hpp"
#include "celestial_body.hpp"
#include <armadillo>
#include <vector>

SolarSystem::SolarSystem() {
  m_kinetic_energy = 0;
  m_potential_energy = 0;
  m_angular_momentum = arma::zeros(3);
}

CelestialBody& SolarSystem::createCelestialBody(arma::vec& pos, arma::vec& vel, double m){
  m_bodies.push_back(CelestialBody(pos,vel,m));
  return m_bodies.back();
}

int SolarSystem::numberOfBodies() const {
  return m_bodies.size();
}

void SolarSystem::calculateForcesAndEnergy() {
  m_kinetic_energy = 0;
  m_potential_energy = 0;
  const double G = 4*M_PI*M_PI;
  m_angular_momentum = arma::zeros(3);

  for (CelestialBody &body: m_bodies) {
    body.force.zeros();
  }

  for (int i = 0; i<numberOfBodies();++i){
    CelestialBody &body1 = m_bodies[i];
    for (int j = i+1; j<numberOfBodies();++j) {
      CelestialBody &body2 = m_bodies[j];

      arma::vec dr_vec = body2.position - body1.position;

      double dr = arma::norm(dr_vec);
      double potential_energy = G*body1.mass*body2.mass/dr;

      arma::vec gravforce = dr_vec*potential_energy/(dr*dr);

      body1.force += gravforce;
      body2.force -= gravforce;

      m_potential_energy += potential_energy;
    }

    double v = arma::norm(body1.velocity);
    m_kinetic_energy += 0.5*body1.mass*v*v;
    m_angular_momentum += arma::cross(body1.position,body1.mass*body1.velocity);
  }
}

std::vector<CelestialBody> &SolarSystem::bodies(){
  return m_bodies;
}

double SolarSystem::potentialEnergy(){
  return m_potential_energy;
}

double SolarSystem::kineticEnergy(){
  return m_kinetic_energy;
}

arma::vec SolarSystem::angularMomentum(){
  return m_angular_momentum;
}

void SolarSystem::initiateDataFile(std::string positions_filename) {
  m_positions_filename = positions_filename;

  m_file.open(m_positions_filename.c_str(), std::ofstream::out);

  int numBods = numberOfBodies();
  m_file << "Number of bodies: " << numBods << ". Number of columns: "
         << 3*numBods << "." << std::endl;

  m_file << "There are three columns per body. The columns are x0 y0 z0 ... x"
         << numBods-1 << " y"
         << numBods-1 << " z"
         << numBods-1 << std::endl;

  m_file << "Rows correspond to time-values t0, t1, ... in descending order."
         << std::endl;
}

void SolarSystem::initiateDataFile(std::string positions_filename, std::string energy_filename) {
  initiateDataFile(positions_filename);
  m_energy_filename = energy_filename;

  m_energy_file.open(m_energy_filename.c_str(),std::ofstream::out);

  m_energy_file << "The five columns in this file are (in order): kinetic energy, potential energy, angular momentum (x,y,z)";
  m_energy_file << std::endl;
}

void SolarSystem::writeToFile() {
  if(!m_file.good()) {
    m_file.open(m_positions_filename.c_str(), std::ofstream::out);
    if(!m_file.good()) {
      std::cout << "Error opening file " << m_positions_filename << ". Aborting!" << std::endl;
      std::terminate();
    }
  }

  for(CelestialBody &body : m_bodies) {
    m_file << body.position(0) << " " << body.position(1) << " " << body.position(2) << " ";
  }
  m_file << std::endl;
}

void SolarSystem::writeEnergyToFile() {
  if(!m_energy_file.good()) {
    m_energy_file.open(m_energy_filename.c_str(), std::ofstream::out);
    if(!m_file.good()) {
      std::cout << "Error opening file " << m_positions_filename << ". Aborting!" << std::endl;
      std::terminate();
    }
  }
  m_energy_file << kineticEnergy() << " " << potentialEnergy() << " " << m_angular_momentum(0) << " " <<
                   m_angular_momentum(1) << " " << m_angular_momentum(2) << std::endl;
}
