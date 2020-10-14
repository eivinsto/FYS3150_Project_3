#include "catch.hpp"
#include "solar_integrator.hpp"
#include "solar_system.hpp"
#include <armadillo>
#include <iostream>

TEST_CASE("Test energy conservation with VelocityVerlet") {
  int numTimesteps = 100000;
  double dt = 0.001;
  arma::vec energy(numTimesteps);

  SolarSystem solarSystem;
  arma::vec suninit = {0,0,0};
  solarSystem.createCelestialBody( suninit, suninit, 1.0 );

  arma::vec earthposinit = {1,0,0};
  arma::vec earthvelinit = {0, 2*M_PI, 0};
  solarSystem.createCelestialBody( earthposinit, earthvelinit, 3e-6 );

  solar_integrator integrator(dt, "VelocityVerlet");

  for(int i=0; i<numTimesteps; i++) {
    integrator.integrateOneStep(solarSystem);
    energy[i] = solarSystem.potentialEnergy();
  }
  REQUIRE(energy[0] == Approx(energy[numTimesteps-1]));
}

TEST_CASE("Test Angular momentum conservation with VelocityVerlet") {
  int numTimesteps = 100000;
  double dt = 0.001;
  arma::vec angmom(numTimesteps);

  SolarSystem solarSystem;
  arma::vec suninit = {0,0,0};
  solarSystem.createCelestialBody( suninit, suninit, 1.0 );

  arma::vec earthposinit = {1,0,0};
  arma::vec earthvelinit = {0, 2*M_PI, 0};
  solarSystem.createCelestialBody( earthposinit, earthvelinit, 3e-6 );

  solar_integrator integrator(dt, "VelocityVerlet");

  for(int i=0; i<numTimesteps; i++) {
    integrator.integrateOneStep(solarSystem);
    arma::vec angmomvec = solarSystem.angularMomentum();
    angmom[i] = angmomvec[0]*angmomvec[0] +
                       angmomvec[1]*angmomvec[1] +
                       angmomvec[2]*angmomvec[2];
  }
  REQUIRE(angmom[0] == Approx(angmom[numTimesteps-1]));
}

TEST_CASE("Test energy conservation with Euler") {
  int numTimesteps = 100000;
  double dt = 0.001;
  arma::vec energy(numTimesteps);

  SolarSystem solarSystem;
  arma::vec suninit = {0,0,0};
  solarSystem.createCelestialBody( suninit, suninit, 1.0 );

  arma::vec earthposinit = {1,0,0};
  arma::vec earthvelinit = {0, 2*M_PI, 0};
  solarSystem.createCelestialBody( earthposinit, earthvelinit, 3e-6 );

  solar_integrator integrator(dt, "Euler");

  for(int i=0; i<numTimesteps; i++) {
    integrator.integrateOneStep(solarSystem);
    energy[i] = solarSystem.potentialEnergy();
  }
  REQUIRE(energy[0] == Approx(energy[numTimesteps-1]));
}

TEST_CASE("Test Angular momentum conservation with Euler") {
  int numTimesteps = 100000;
  double dt = 0.001;
  arma::vec angmom(numTimesteps);

  SolarSystem solarSystem;
  arma::vec suninit = {0,0,0};
  solarSystem.createCelestialBody( suninit, suninit, 1.0 );

  arma::vec earthposinit = {1,0,0};
  arma::vec earthvelinit = {0, 2*M_PI, 0};
  solarSystem.createCelestialBody( earthposinit, earthvelinit, 3e-6 );

  solar_integrator integrator(dt, "Euler");

  for(int i=0; i<numTimesteps; i++) {
    integrator.integrateOneStep(solarSystem);
    arma::vec angmomvec = solarSystem.angularMomentum();
    angmom[i] = angmomvec[0]*angmomvec[0] +
                       angmomvec[1]*angmomvec[1] +
                       angmomvec[2]*angmomvec[2];

  }
  REQUIRE(angmom[0] == Approx(angmom[numTimesteps-1]));
}
