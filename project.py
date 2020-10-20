"""
Python script to interface with project code.
"""
from subprocess import run
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import os
import itertools


# retriveing working directories:
rootdir = os.getcwd()
src = rootdir + "/src"


def build_cpp():
    """Function for building c++ program."""
    run(["make", "all"], cwd=src)


def test_cpp():
    """Function for running unit-tests."""
    run(["make", "test"], cwd=src)
    run("./test_main.exe", cwd=src)


class SolarSystem:
    """SolarSystemFiles is a class for reading data from SolarSystem C++
    program, and plotting orbits and energy/momentum."""

    def __init__(self, numTimesteps, dt, write_limit, integration_method,
                 init_file, posfile, momenfile, bodynames):
        """Class constructor.
        Args:
            numTimesteps: integer - number of time steps to take.
            dt: float (64bit) - size of time step.
            write_limit: integer - Number of time steps between each
                                   write to file. Example: write_limit = 10
                                   means data is written every 10 time steps.
            integration_method: string - name of the method of integration to
                                use 'Euler' for forward Euler, and
                                'VelocityVerlet' for Velocity-Verlet.
            init_file: string - name of file to load inital conditions from.
            posfile: string - Name of datafile containing plannet positions.
            momenfile: string - Name of datafile containing energy and angular
                                momentum.
            bodynames: list of strings - List of names of each celestial body.
        """
        self.isgenerated = False

        self.numTimesteps = numTimesteps
        self.dt = dt
        self.integration_method = integration_method
        # adding path to filenames:
        self.init_file = rootdir + "/data/" + init_file
        self.posfile = rootdir + "/data/" + posfile
        self.momenfile = rootdir + "/data/" + momenfile

        self.bodynames = bodynames

        # ensuring that maximum 1000 lines of data are read from each file,
        # this is to limit memory use.
        if (self.numTimesteps//write_limit > 1000):
            self.everyNlines = self.numTimesteps // (write_limit*1000)
        else:
            self.everyNlines = 1

        # generating array of Simulated time:
        self.times = np.linspace(0, numTimesteps,
                                 numTimesteps//(self.everyNlines*write_limit))

    def generateSystem(self):
        """Method for calling C++ program and generating data."""
        build_cpp()  # building c++ program if necessary.

        # passing arguments to c++ program and running simulation:
        run(
            [
                "./main.exe",
                f"{numTimesteps}",
                f"{dt}",
                f"{write_limit}",
                self.integration_method,
                self.init_file,
                self.posfile,
                self.momenfile,
            ],
            cwd=src
        )

        self.readData()  # reading data from files
        print("Done generating data")
        self.isgenerated = True  # setting flag that data has been generated.

    def readData(self):
        """Method for reading data from files."""

        # reading header and setting number of bodies:
        with open(self.posfile) as infile:
            header1 = infile.readline().split()
            word = header1[3].rstrip('.')
            self.numBods = int(word)

            # reading position data from file:
            self.bodyPos = np.genfromtxt(itertools.islice(
                infile, 2, self.numTimesteps+2, self.everyNlines))
            # self.bodyPos is an Nx3n-array where each row contains the
            # positions of the bodies in a moment in time.
            # the positions are stored as coulumns in the patern
            # x_0 y_0 z_0 x_1 y_1 z_1 ... x_n-1 y_n-1 z_n-1
            # where n is the number of bodies.

        # reading angular momentum from file:
        with open(self.momenfile, "r") as infile:
            self.angmom = np.genfromtxt(itertools.islice(
                infile, 1, self.numTimesteps+1, self.everyNlines),
                usecols=[2, 3, 4])
            # self.angmom is an Nx3-array containing the components of the
            # angular momentum of the System. Each row is a single moment
            # in time, and the coulumns are x y z.

        # reading Kinetic and Potential energy from file:
        with open(self.momenfile, "r") as infile:
            self.energy = np.genfromtxt(itertools.islice(
                infile, 1, self.numTimesteps+1, self.everyNlines),
                usecols=[0, 1])
            # self.energy is an Nx2-array containing the Kinetic energy for
            # each time step as the first column and potential energy
            # each time step as the second column.

    def orbit2D(self, number_of_bodies=None):
        """Method for creating 2D plot of orbits,
           by plotting the x and y coordinates of the bodies' positions.
           Args:
               number_of_bodies: int - number of bodies to plot,
                                       counting from zero"""
        if not self.isgenerated:  # generating data if necessary.
            self.generateSystem()
        if number_of_bodies is None:
            number_of_bodies = self.numBods

        plt.figure()  # creates figure

        # running through celestial bodies:
        for i in range(number_of_bodies):
            if i == 0:
                plottype = "."  # marking sun with dotted orbit.
            else:
                plottype = "-"

            # plotting orbits relative to center of mass:
            plt.plot(self.bodyPos[:, 3*i],
                     self.bodyPos[:, 3*i + 1],
                     plottype,
                     label=self.bodynames[i])

        # adding titles and legend:
        plt.title(
            f"n = {self.numBods}, N = {self.numTimesteps:.1e}, dt = {self.dt},"
            + " " + integration_method +
            f", Simulated time = {self.dt*self.numTimesteps} years"
        )
        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        plt.legend()
        plt.grid()
        plt.axis('equal')

    def orbit3D(self, number_of_bodies=None):
        """Create 3D plot of orbits.
        Args:
            number_of_bodies: int - number of bodies to plot,
                                    counting from zero"""
        if not self.isgenerated:  # generating data if necessary.
            self.generateSystem()
        if number_of_bodies is None:
            number_of_bodies = self.numBods

        fig = plt.figure()  # creates figure
        ax = Axes3D(fig)  # create 3D subplot

        # running through celestial bodies:
        for i in range(number_of_bodies):
            if i == 0:
                plottype = "."  # marking sun with dotted orbit.
            else:
                plottype = "-"

            # plotting orbit of current body:
            ax.plot3D(self.bodyPos[:, 3*i],
                      self.bodyPos[:, 3*i + 1],
                      self.bodyPos[:, 3*i + 2],
                      plottype,
                      label=self.bodynames[i])

        # adding labels:
        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        ax.set_zlabel("z [AU]")
        ax.legend()

    def plotEnergy(self):
        """Plotting kinetic, potential and total energy for system as
        function of time steps. Values are scaled by maximum magnitude of
        total energy."""
        if not self.isgenerated:  # generating data if necessary.
            self.generateSystem()

        totenergy = self.energy[:, 0] + self.energy[:, 1]
        emax = np.max(np.abs(totenergy))
        totenergy_mean = np.mean(totenergy/emax)
        totenergy_std = np.std(totenergy/emax)

        plt.figure()  # creating figure

        # plotting kinetic, potential and total energy
        plt.plot(self.times, self.energy[:, 0]/emax,
                 label="Kinetic energy")
        plt.plot(self.times, self.energy[:, 1]/emax,
                 label="Potential energy")
        plt.plot(self.times, totenergy/emax, '--',
                 label=r"Total energy. E_mean = " +
                 f"{totenergy_mean:.2e}, E_std = {totenergy_std:.2e}")

        # adding titles and lables:
        plt.title(
            f"n = {self.numBods}, N = {self.numTimesteps:.1e}, dt = {self.dt},"
            + " " + integration_method +
            f", Simulated time = {self.dt*self.numTimesteps} years"
        )
        plt.xlabel("N [number of time steps]")
        plt.ylabel("Energy $[E_{tot}/max(|E_{tot}|)]$")
        plt.legend()
        plt.grid()

    def plotAngMomMagnitude(self):
        if not self.isgenerated:  # generating data if necessary.
            self.generateSystem()

        # calculating magnitude of angular momentum:
        angmommag = np.sqrt(self.angmom[:, 0]**2 +
                            self.angmom[:, 1]**2 +
                            self.angmom[:, 2]**2)
        Lmax = np.max(np.abs(angmommag))
        Lmean = np.mean(angmommag/Lmax)
        Lstd = np.std(angmommag/Lmax)

        plt.figure()  # creating figure

        # plotting magnitude of angular momentum scaled by maximum magnitude:
        plt.plot(self.times, angmommag/Lmax,
                 label="Angular momentum magnitude. L_mean = " +
                 f"{Lmean:.2e}, L_std = {Lstd:.2e}")

        # adding titles and legend:
        plt.title(
            f"n = {self.numBods}, N = {self.numTimesteps:.1e}, dt = {self.dt},"
            + " " + integration_method +
            f", Simulated time = {self.dt*self.numTimesteps} years"
        )
        plt.xlabel("N [number of time steps]")
        plt.ylabel("Angular momentum $[L_{tot}/max(|L_{tot}|)]$")
        plt.legend()
        plt.grid()


# name of bodies used in project:
bodynames = ["Sun", "Mercury", "Venus", "Earth", "Mars",
             "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

# asking for imput:
print("""Write se for Sun-Earth simulation,
sej for Sun-Earth-Jupiter,
sm for Sun-Mercury,
and ss for entire Solar System.
Write test to run unit-tests,
or b for benchmark.""")
runflag = input("Choose run: ")

if runflag != "test" or runflag != "b":
    # asking for imput parameters from user:
    numTimesteps = int(eval(input("Number of time steps N = ")))
    dt = float(eval(input("Size of time step dt = ")))
    limit_write = input(
        "Only write the data from 1000 evenly spaced time steps? y/n: "
    )

if runflag == "se":  # initial data for sun_earth run:
    init_file = "earth-sun-init.txt"
    bodynames = [bodynames[0], bodynames[3]]

elif runflag == "sej":  # initial data for sun_earth_jupiter run:
    init_file = "sun-earth-jupiter-2020-Oct-19-00:00:00.init"
    bodynames = [bodynames[0], bodynames[3], bodynames[5]]

elif runflag == "sm":  # initial data for sun_mercury run:
    init_file = "sun-mercury.init"
    bodynames = [bodynames[0], bodynames[1]]

elif runflag == "ss":  # initial data for entire SolarSystem run:
    init_file = "sun-and-friends-2020-Oct-19-00:00:00.init"

if runflag != "test" or runflag != "b":  # setting up run:
    if limit_write == "y":
        write_limit = numTimesteps//1000
    else:
        write_limit = 1

    integration_method = input("""Choose integration method:
    Write fe for forward Euler,
    or vv for Velocity-Verlet: """)

    # setting up algorithm and output filenames:
    if integration_method == "vv":
        integration_method = "VelocityVerlet"
    elif integration_method == "fe":
        integration_method = "Euler"

    posfile = runflag + "_" + integration_method + "_" + "positions.xyz"
    momenfile = runflag + "_" + integration_method + "_" + "energies.dat"

    # initalising instance of SolarSystem class with parameters:
    sun_earth = SolarSystem(
        numTimesteps,
        dt,
        write_limit,
        integration_method,
        init_file,
        posfile,
        momenfile,
        bodynames
    )

    # calling the various plot functions:
    sun_earth.orbit3D()
    sun_earth.orbit2D()
    sun_earth.plotEnergy()
    sun_earth.plotAngMomMagnitude()

    plt.show()

elif runflag == "test":
    test_cpp()
