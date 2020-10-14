#define _USE_MATH_DEFINES

#include <cmath>
#include "solar_system.hpp"
#include "celestial_body.hpp"
#include <armadillo>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>


/**
* Constructor that takes only the exponent to be used in the force as argument.
* SolarSystem object is created empty in this case.
*
* @beta Exponent in the force
*/
SolarSystem::SolarSystem(double beta) {
  m_kinetic_energy = 0;
  m_potential_energy = 0;
  m_angular_momentum = arma::zeros(3);
  m_beta = beta;
}

/**
* Basic constructor that takes no argument. Creates an empty SolarSystem object
* with the exponent in the force set to 2.
*/
SolarSystem::SolarSystem() : SolarSystem(2.0) {}

/**
* Constructor that takes an input file and reads initial conditions from said file.
* A parameter beta also specifies the exponent in the force. Initial conditions of one object
* should be on a single line: first the three position coordinates (x,y,z), then
* three velocity coordinates (vx,vy,vz) and the last element should be the mass.
* Elements should be separated by whitespace.
* Position should be in units AU, velocity in AU y^-1 and mass in units of the
* mass of the sun.
*
* @input_filename Name of file to read initial conditions from.
* @beta Exponent in the force
*/
SolarSystem::SolarSystem(std::string input_filename, double beta) : SolarSystem(beta) {
  std::ifstream input_file(input_filename); // Get input stream
  // Define variables to be used in iteration
  std::string line;
  std::string word;
  arma::vec x = arma::zeros(3);
  arma::vec v = arma::zeros(3);
  double mass = 0;

  // Iterate over the lines in the file
  while ( std::getline(input_file, line)) {
    // Clear variables
    word = "";
    x = arma::zeros(3);
    v = arma::zeros(3);
    m = 0;

    // Redefine line as a new stringstream variable ssline so that getline() can be used to split around whitespaces.
    std::stringstream ssline(line);

    // Read positions
    for (int i = 0; i<3; ++i){
      std::getline(ssline,word,' ');
      x(i) = std::stod(word);
    }

    // Read velocities
    for (int i = 0; i<3; ++i){
      std::getline(ssline,word,' ');
      v(i) = std::stod(word);
    }

    // Read mass
    std::getline(ssline,word,' ');
    mass = std::stod(word);

    // Generate CelestialBody object and add to current SolarSystem object.
    createCelestialBody(x,v,mass);
  }
}

/**
* Constructor that only takes an input file name and calls
* SolarSystem(std::string input_filename, double beta) with said file name and
* beta = 2.
*
* @input_filename Name of file to read initial conditions from.
*/
SolarSystem::SolarSystem(std::string input_filename) : SolarSystem(input_filename, 2.0){}

/**
* Add a new CelestialBody object to the SolarSystem object.
*
* @pos Vector containing position of the new body.
* @vel Vector containing position of the new body.
* @m Double containing mass of the new body.
*/
CelestialBody& SolarSystem::createCelestialBody(arma::vec& pos, arma::vec& vel, double m){
  m_bodies.push_back(CelestialBody(pos,vel,m));
  return m_bodies.back();
}

/**
* Returns number of CelestialBody objects currently in the SolarSystem object.
*/
int SolarSystem::numberOfBodies() const {
  return m_bodies.size();
}

/**
* Function that calculates kinetic and potential energy, angular momentum,
* and the force that the bodies exert on each other.
*/
void SolarSystem::calculateForcesAndEnergy() {
  m_kinetic_energy = 0;
  m_potential_energy = 0;
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
      double potential_energy = m_G*body1.mass*body2.mass/dr;

      arma::vec gravforce = dr_vec*potential_energy/std::pow(dr,m_beta);

      body1.force += gravforce;
      body2.force -= gravforce;

      m_potential_energy += potential_energy;
    }

    double v = arma::norm(body1.velocity);
    m_kinetic_energy += 0.5*body1.mass*v*v;
    m_angular_momentum += arma::cross(body1.position,body1.mass*body1.velocity);
  }
}

void SolarSystem::calculateForcesAndEnergyWithRelativisticCorrection() {
  m_kinetic_energy = 0;
  m_potential_energy = 0;
  m_angular_momentum = arma::zeros(3);


  for (CelestialBody &body: m_bodies) {
    body.force.zeros();
  }

  for (int i = 0; i<numberOfBodies();++i){
    CelestialBody &body1 = m_bodies[i];
    arma::vec l_vec = arma::cross(body1.position,body1.mass*body1.velocity);
    double l = arma::norm(l_vec);
    double l2 = l*l;


    for (int j = i+1; j<numberOfBodies();++j) {
      CelestialBody &body2 = m_bodies[j];

      arma::vec dr_vec = body2.position - body1.position;

      double dr = arma::norm(dr_vec);
      double dr2 = dr*dr;
      double potential_energy = m_G*body1.mass*body2.mass/dr;


      arma::vec gravforce = dr_vec*(potential_energy/dr2) * (1 + m_rel_constant*l2/dr2 );

      body1.force += gravforce;
      body2.force -= gravforce;

      m_potential_energy += potential_energy;
    }

    double v = arma::norm(body1.velocity);
    m_kinetic_energy += 0.5*body1.mass*v*v;
    m_angular_momentum += l_vec;
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
