#include "celestial_body.hpp"
#include <armadillo>

CelestialBody::CelestialBody(arma::vec& pos, arma::vec& vel, double m) {
  position = pos;
  velocity = vel;
  mass = m;
  force = arma::zeros(3);
}
