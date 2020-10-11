#ifndef SOLAR_SYSTEM_HPP
#define SOLAR_SYSTEM_HPP

#include "celestial_body.hpp"
#include <vector>

class SolarSystem {
public:
  SolarSysem();
  createCelestialBody(vec& pos, vec& vel, double m);

private:
  std::vector<CelestialBody> bodies;
};

#endif
