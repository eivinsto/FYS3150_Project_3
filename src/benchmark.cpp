#include "catch.hpp"
#include "solar_integrator.hpp"
#include "solar_system.hpp"
#include <time.h>
#include <string>
#include <armadillo>
#include <iostream>

int main() {
  int numTimesteps = 1000000;  // number of time steps
  double dt = 0.000248;  // size of time step
  int numRuns = 50;  // number of simulations to run
  arma::mat benchmark_times(numRuns, 2);  // array for saving timings

  // declaring clock_t objects for timing
  clock_t start, finish;

  for (int j = 0; j < numRuns; j++) {
    std::cout << "\tPerforming benchmark of forward Euler: " << 100*(j+1)/numRuns << "%\r";
    std::cout.flush();

    // creating simulation:
    std::string init_file = "../data/sun-and-friends-2020-Oct-19-00:00:00.init";  // File containing initial conditions

    // Reading initial state from file
    SolarSystem solarSystem(init_file);
    solar_integrator integrator(dt, "Euler");

    // timing integration
    start = clock();
    for(int i=0; i<numTimesteps; i++) {
      integrator.integrateOneStep(solarSystem);
    }
    finish = clock();

    // saving result to matrix
    benchmark_times(j,0) = double(finish - start)/CLOCKS_PER_SEC;
  }
  std::cout << std::endl;

  for (int j = 0; j < numRuns; j++) {
    std::cout << "\tPerforming benchmark of VelocityVerlet: " << 100*(j+1)/numRuns << "%\r";
    std::cout.flush();

    // creating simulation:
    std::string init_file = "../data/sun-and-friends-2020-Oct-19-00:00:00.init";  // File containing initial conditions

    // Reading initial state from file
    SolarSystem solarSystem(init_file);
    solar_integrator integrator(dt, "VelocityVerlet");

    // timing integration
    start = clock();
    for(int i=0; i<numTimesteps; i++) {
      integrator.integrateOneStep(solarSystem);
    }
    finish = clock();

    // saving result to matrix
    benchmark_times(j,1) = double(finish - start)/CLOCKS_PER_SEC;
  }
  std::cout << std::endl;

  // saving result to file
  benchmark_times.save("../data/benchmarkdata.dat", arma::arma_ascii);
  return 0;
}
