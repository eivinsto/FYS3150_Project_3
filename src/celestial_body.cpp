#include "celestial_body.hpp"

CelestialBody::CelestialBody(arma::vec& pos, arma::vec& vel, double m) {
  position = pos;
  velocity = vel;
  mass = m;
}
