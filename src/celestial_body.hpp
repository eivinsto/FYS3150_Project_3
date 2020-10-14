#ifndef CELESTIAL_BODY_HPP
#define CELESTIAL_BODY_HPP

#include <armadillo>

class CelestialBody {
public:
  // Internally stored variables
  arma::vec position;
  arma::vec velocity;
  double mass;
  arma::vec force;

  // Constructor
  CelestialBody(arma::vec& pos, arma::vec& vel, double m);
};

#endif
