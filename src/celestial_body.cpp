#include "celestial_body.hpp"

CelestialBody::CelestialBody(vec& pos, vec& vel, double m) {
  position = pos;
  velocity = vel;
  mass = m;
}
