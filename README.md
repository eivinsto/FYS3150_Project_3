# FYS3150_Project_3
All code for this report was written in python 3.8, and C++.
To generate the data used in the report, run the python-script "project.py" in the root directory of this repo.

## Usage:
When executed, the "project.py" script will ask you which part of the report to run:
*   se   - Run Sun-Earth simulation. There will be additional options for the different simulations using the Sun-Earth system.
*   sej  - Run simulations using Sun-Earth-Jupiter system.
*   sm   - Run simulations of perihelion precession using Sun-Mercury system.
*   ss   - Run simulation of Sun with all planets and Pluto.
*   test - Run unit tests for the solar\_integrator algorithms.
*   b    - Run benchmark of solar\_integrator algorithms simulating the Sun with all planets and Pluto.

For all simulations, except for test and benchmark, you will be asked to input the number of time steps, and size of time step in years. You will also be asked if data should be saved to file on every time step, or be limited to 10000 evenly spaced time steps. Choosing to limit the number of writes speeds up integration:
```console
$ python project.py
Write se for Sun-Earth simulation,
sej for Sun-Earth-Jupiter,
sm for Sun-Mercury,
and ss for entire Solar System.
Write test to run unit-tests,
or b for benchmark.
Choose run: se

Number of time steps N = 1e7
Size of time step dt = 1e-7

Only write the data from 10000 evenly spaced time steps? y/n: y
```
##### Sun-Earth:
For the Sun-Earth system you will be given the choice to run an escape velocity simulation where you can input the initial, tangential velocity of the Earth. The input method supports numbers, as well as python expressions that evaluate to numbers. Thus you can give the exact value:

```console
Run escape velocity test? y/n: y
Enter initial escape velocity: 2*np.pi*np.sqrt(2)
```
Declining the escape velocity simulation, prompts the simulations with varying beta:
```console
Run escape velocity test? y/n: n
Run sun-earth simulation with varying beta? y/n: y
Use elliptical orbit? y/n: y
```
The "Use elliptical orbit?" option lets you choose if the Earth should start with initial conditions that give circular orbit (n) or eccentric orbit (y), for the simulation with varying beta.

##### Sun with all planets and Pluto:
For the simulation with the Sun and all planets + Pluto the option will be given to only plot the orbits of the Sun and the 5 innermost planets, for easier viewing of their orbits in the 3D plot.

Finally, for any run the script will prompt you to choose which algorithm to integrate with:
*   fe - Forward Euler.
*   vv - Velocity-Verlet.

```console
Choose integration method:
    Write fe for forward Euler,
    or vv for Velocity-Verlet:
    fe
```

##### Example run of Sun-Mercury simulation:
```console
$ python project.py
Write se for Sun-Earth simulation,
sej for Sun-Earth-Jupiter,
sm for Sun-Mercury,
and ss for entire Solar System.
Write test to run unit-tests,
or b for benchmark.
Choose run: sm

Number of time steps N = 1e8
Size of time step dt = 1e-6

Choose integration method:
    Write fe for forward Euler,
    or vv for Velocity-Verlet:
    vv
g++ -Wall -Wextra -O2 -march=native -g main.cpp -c
g++ -Wall -Wextra -O2 -march=native -g main.o solar_integrator.o solar_system.o celestial_body.o -o main.exe -larmadillo
Current beta = 2

Done generating data
make: Nothing to be done for 'all'.
Current beta = 2

Done generating data
ThetaP without correction = 0.681256''
ThetaP with correction = 43.0864''
```
