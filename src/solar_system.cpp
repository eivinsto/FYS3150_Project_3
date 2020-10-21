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
    mass = 0;

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
  moveToCOFMFrame();
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
  // Clear variables
  m_kinetic_energy = 0;
  m_potential_energy = 0;
  m_angular_momentum = arma::zeros(3);
  for (CelestialBody &body: m_bodies) {
    body.force.zeros();
  }

  // Iterate over every body (body1)
  for (int i = 0; i<numberOfBodies();++i){
    CelestialBody &body1 = m_bodies[i];
    // Second iteration over every body after body1 (body2)
    for (int j = i+1; j<numberOfBodies();++j) {
      CelestialBody &body2 = m_bodies[j];

      // Distance vector
      arma::vec dr_vec = body2.position - body1.position;

      // Length of distance vector
      double dr = arma::norm(dr_vec);

      // Potential energy between body1 and body2
      double potential_energy = m_G*body1.mass*body2.mass/std::pow(dr,m_beta-1);

      // Force vector between the bodies (sign adjusted when adding to the bodies)
      arma::vec gravforce = dr_vec*potential_energy/(dr*dr);
      body1.force += gravforce;
      body2.force -= gravforce;

      // Adding potential energy to total
      m_potential_energy -= potential_energy;
    }

    // Adding kinetic energy of body1 to total
    double v = arma::norm(body1.velocity);
    m_kinetic_energy += 0.5*body1.mass*v*v;

    // Adding angular momentum of body1 to total angular momentum
    m_angular_momentum += arma::cross(body1.position,body1.mass*body1.velocity);
  }
}

/**
* Function that calculates kinetic and potential energy, angular momentum,
* and the force that the bodies exert on each other with a relativistic correction.
*/
void SolarSystem::calculateForcesWithRelativisticCorrection() {
  // Clear variables
  for (CelestialBody &body: m_bodies) {
    body.force.zeros();
  }

  // Load Sun (body1)
  CelestialBody &body1 = m_bodies[0];

  // Load Mercury (body2)
  CelestialBody &body2 = m_bodies[1];

  // Calculate angular momentum vector, size and squared size
  arma::vec l_vec = arma::cross(body2.position, body2.velocity);
  double l = arma::norm(l_vec);
  double l2 = l*l;


  // Distance vector
  arma::vec dr_vec = body2.position;

  // Length of distance vector
  double dr = arma::norm(dr_vec);
  double dr2 = dr*dr;

  // Potential energy between body1 and body2
  double potential_energy = m_G*body1.mass*body2.mass/dr;

  // Calculate force vector with relativistic correction factor (sign adjusted when adding to the bodies)
  arma::vec gravforce = dr_vec*(potential_energy/dr2) * (1 + m_rel_constant*l2/dr2 );
  body2.force -= gravforce;
}

/**
* Returns a vector containing the CelestialBody objects currently in the
* SolarSystem obect.
*/
std::vector<CelestialBody> &SolarSystem::bodies(){
  return m_bodies;
}

/**
* Returns total potential energy. If calculateForcesAndEnergy() has not been run
* this returns 0.
*/
double SolarSystem::potentialEnergy(){
  return m_potential_energy;
}

/**
* Returns total kinetic energy. If calculateForcesAndEnergy() has not been run
* this returns 0.
*/
double SolarSystem::kineticEnergy(){
  return m_kinetic_energy;
}

/**
* Returns total angular momentum. If calculateForcesAndEnergy() has not been run
* this returns a 0-vector.
*/
arma::vec SolarSystem::angularMomentum(){
  return m_angular_momentum;
}

/**
* This function initiates a data file to store the positions of the objects.
* It should be followed up by using writeToFile() in every time step to actually
* write data in the file.
*
* @positions_filename Name of file to store position data in.
*/
void SolarSystem::initiateDataFile(std::string positions_filename) {
  // Save filename internally
  m_positions_filename = positions_filename;

  // Open file
  m_file.open(m_positions_filename.c_str(), std::ofstream::out);

  // Write header in file describing contents
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

/**
* This function initiates a data file to store the positions of the objects, and
* another file in which to store total kinetic and potential energy, and angular
* momentum.
* It should be followed up by using writeEnergyToFile() in every time step to
* actually write data in the file.
*
* @positions_filename Name of file to store position data in.
* @energy_filename Name of file to store energies and angular momentum in.
*/
void SolarSystem::initiateDataFile(std::string positions_filename, std::string energy_filename) {
  // Inititate datafile for the positions
  initiateDataFile(positions_filename);

  // Store filename for energies and angular momentum internally
  m_energy_filename = energy_filename;

  // Open file
  m_energy_file.open(m_energy_filename.c_str(),std::ofstream::out);

  // Write header describing contents of file.
  m_energy_file << "The five columns in this file are (in order): kinetic energy, potential energy, angular momentum (x,y,z)";
  m_energy_file << std::endl;
}

/**
* Function that writes positions to file in current timestep. The file must be
* initialized using initiateDataFile() before this function can be used.
*/
void SolarSystem::writeToFile() {
  // Check whether file is still working as it should
  if(!m_file.good()) {
    m_file.open(m_positions_filename.c_str(), std::ofstream::out);
    if(!m_file.good()) {
      std::cout << "Error opening file " << m_positions_filename << ". Aborting!" << std::endl;
      std::terminate();
    }
  }

  // Write positions to file
  for(CelestialBody &body : m_bodies) {
    m_file << body.position(0) << " " << body.position(1) << " " << body.position(2) << " ";
  }

  // Change line
  m_file << std::endl;
}

/**
* Function that writes energies and angular momentum in current timestep
* to file. The file must be initialized using initiateDataFile() before this
* function can be used.
*/
void SolarSystem::writeEnergyToFile() {
  // Check whether file is still working as it should
  if(!m_energy_file.good()) {
    m_energy_file.open(m_energy_filename.c_str(), std::ofstream::out);
    if(!m_file.good()) {
      std::cout << "Error opening file " << m_positions_filename << ". Aborting!" << std::endl;
      std::terminate();
    }
  }

  // Write energies and angular momentum to file
  m_energy_file << kineticEnergy() << " " << potentialEnergy() << " " << m_angular_momentum(0) << " " <<
                   m_angular_momentum(1) << " " << m_angular_momentum(2) << std::endl;
}

/**
* Function that moves the system to the center of mass frame.
*/
void SolarSystem::moveToCOFMFrame(){
  // Sum weighted positions, momentum and mass in r, v, and tot_mass respectively
  arma::vec r = arma::zeros(3);
  arma::vec v = arma::zeros(3);
  double tot_mass = 0;
  for (int i = 0; i<numberOfBodies(); ++i){
    CelestialBody& body = m_bodies[i];
    r += body.mass*body.position;
    v += body.mass*body.velocity;
    tot_mass += body.mass;
  }

  // Get center of mass position and velocity (overwrite r and v)
  r = r/tot_mass;
  v = v/tot_mass;

  // Subtract center of mass position and velocity from all positions and velocity
  for (int i = 0; i<numberOfBodies(); ++i){
    CelestialBody& body = m_bodies[i];
    body.position -= r;
    body.velocity -= v;
  }
}
