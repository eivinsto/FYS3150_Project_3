#include "celestial_body.hpp"
#include <armadillo>

/**
* Constructor that takes a position vector and a velocity vector and a double mass.
* These are stored internally in the class.
*/
CelestialBody::CelestialBody(arma::vec& pos, arma::vec& vel, double m) {
  position = pos;
  velocity = vel;
  mass = m;

  // Sets force to zero so it is initialized.
  force = arma::zeros(3);
}
