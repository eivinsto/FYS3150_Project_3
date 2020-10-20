#include "catch.hpp"
#include "solar_integrator.hpp"
#include "solar_system.hpp"
#include <time.h>
#include <armadillo>
#include <iostream>

int main() {
  int numTimesteps = 10000;
  double dt = 0.0001;
  int numRuns = 5000;
  arma::mat benchmark_times(numRuns, 2);

  // declaring clock_t objects for timing
  clock_t start, finish;

  for (int j = 0; j < numRuns; j++) {
    std::cout << "\tPerforming benchmark of forward Euler: " << 100*(j+1)/numTimesteps << "%\r";
    std::cout.flush();

    SolarSystem solarSystem;
    arma::vec suninit = {0,0,0};
    solarSystem.createCelestialBody( suninit, suninit, 1.0 );

    arma::vec earthposinit = {1,0,0};
    arma::vec earthvelinit = {0, 2*M_PI, 0};
    solarSystem.createCelestialBody( earthposinit, earthvelinit, 3e-6 );

    solar_integrator integrator(dt, "Euler");

    start = clock();
    for(int i=0; i<numTimesteps; i++) {
      integrator.integrateOneStep(solarSystem);
    }
    finish = clock();

    // saving result to matrix
    benchmark_times(j,0) = double(finish - start)/CLOCKS_PER_SEC;
  }

  for (int j = 0; j < numRuns; j++) {
    std::cout << "\tPerforming benchmark of VelocityVerlet: " << 100*(j+1)/numTimesteps << "%\r";
    std::cout.flush();

    SolarSystem solarSystem;
    arma::vec suninit = {0,0,0};
    solarSystem.createCelestialBody( suninit, suninit, 1.0 );

    arma::vec earthposinit = {1,0,0};
    arma::vec earthvelinit = {0, 2*M_PI, 0};
    solarSystem.createCelestialBody( earthposinit, earthvelinit, 3e-6 );

    solar_integrator integrator(dt, "VelocityVerlet");

    start = clock();
    for(int i=0; i<numTimesteps; i++) {
      integrator.integrateOneStep(solarSystem);
    }
    finish = clock();

    // saving result to matrix
    benchmark_times(j,1) = double(finish - start)/CLOCKS_PER_SEC;
  }

  benchmark_times.save("../data/benchmarkdata.dat", arma::arma_ascii);
  return 0;
}
