#ifndef CELESTIAL_BODY_HPP
#define CELESTIAL_BODY_HPP

#include <armadillo>

class CelestialBody {
public:
  arma::vec position;
  arma::vec velocity;
  double mass;
  arma::vec force;

  CelestialBody(arma::vec& pos, arma::vec& vel, double m);
};

#endif
