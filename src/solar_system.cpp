#define _USE_MATH_DEFINES

#include <cmath>
#include "solar_system.hpp"
#include "celestial_body.hpp"
#include <armadillo>
#include <vector>

SolarSystem::SolarSystem(){
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

      body1.force -= gravforce;
      body2.force += gravforce;

      m_potential_energy += potential_energy;
    }

    double v = arma::norm(body1.velocity);
    m_kinetic_energy += 0.5*body1.mass*v*v;
    m_angular_momentum += arma::cross(body1.position,body1.mass*body1.velocity);
  }
}

std::vector<CelestialBody> SolarSystem::bodies(){
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

void SolarSystem::writeToFile(std::string filename) {
    if(!m_file.good()) {
        m_file.open(filename.c_str(), std::ofstream::out);
        if(!m_file.good()) {
            std::cout << "Error opening file " << filename << ". Aborting!" << std::endl;
            std::terminate();
        }
    }

    for(CelestialBody &body : m_bodies) {
        m_file << "1 " << body.position(0) << " " << body.position(1) << " " << body.position(2) << " ";
    }
    m_file << std::endl;
}
