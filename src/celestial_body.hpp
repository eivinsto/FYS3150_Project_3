#ifndef CELESTIAL_BODY_HPP
#define CELESTIAL_BODY_HPP

#include <armadillo>

class CelestialBody {
public:
  vec position(3);
  vec velocity(3);
  double mass;
  vec force(3);

  CelestialBody(vec& pos, vec& vel, double m);
};

#endif
