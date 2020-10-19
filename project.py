"""
Python script to interface with project code.
"""
from subprocess import run
import numpy as np
import matplotlib.pyplot as plt
import os
import itertools
from mpl_toolkits.mplot3d import Axes3D

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


def clean(files="dat"):
    """Function for cleaning datafiles in src directory."""
    if files == "dat":
        run(["make", "cleandat"], cwd=src)
    if files == "all":
        run(["make", "clean"], cwd=src)


class SolarSystem:
    """SolarSystemFiles is a class for reading data from SolarSystem C++
    program, and plotting orbits and energy/momentum."""

    def __init__(self, numTimesteps, dt, write_limit, integration_method,
                 init_file, posfile, momenfile, bodynames):
        """Class constructor.
        Args:
            posfile: string - Name of datafile containing plannet positions.
            momenfile: string - Name of datafile containing energy and angular
                                momentum.
            bodynames: list of strings - List of names of each celestial body.
            numTimesteps: integer - number of time steps to take.
            dt: float (64bit) - size of time step.
        """
        self.isgenerated = False

        self.numTimesteps = numTimesteps
        self.dt = dt
        self.integration_method = integration_method
        self.init_file = rootdir + "/data/" + init_file
        self.posfile = rootdir + "/data/" + posfile
        self.momenfile = rootdir + "/data/" + momenfile
        self.bodynames = bodynames

        if (self.numTimesteps > 1000):
            self.everyNlines = self.numTimesteps//1000
        else:
            self.everyNlines = 1

        self.times = np.linspace(0,
                                 numTimesteps*dt,
                                 numTimesteps//self.everyNlines)

    def generateSystem(self):
        build_cpp()
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

        self.readData()
        print("Done generating data")
        self.isgenerated = True

    def readData(self):
        """Read header from datafile and set number of bodies."""

        with open(self.posfile) as infile:
            header1 = infile.readline().split()
            word = header1[3].rstrip('.')
            self.numBods = int(word)

            self.bodyPos = np.genfromtxt(itertools.islice(
                infile, 2, self.numTimesteps+2, self.everyNlines))

        with open(self.momenfile, "r") as infile:
            self.angmom = np.genfromtxt(itertools.islice(
                infile, 1, self.numTimesteps+1, self.everyNlines),
                usecols=[2, 3, 4])

        with open(self.momenfile, "r") as infile:
            self.energy = np.genfromtxt(itertools.islice(
                infile, 1, self.numTimesteps+1, self.everyNlines),
                usecols=[0, 1])

        print(self.energy.shape, self.angmom.shape,
              self.bodyPos.shape, self.numTimesteps, self.everyNlines)

    def orbit2D(self):
        """Create 2D plot of orbits."""
        if not self.isgenerated:
            self.generateSystem()

        plt.figure()  # creates figure

        # running through celestial bodies:
        for i in range(self.numBods):
            if i == 0:
                plottype = "."
            else:
                plottype = "-"

            plt.plot(self.bodyPos[:, 3*i],
                     self.bodyPos[:, 3*i + 1],
                     plottype,
                     label=self.bodynames[i])

        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        plt.legend()
        plt.grid()
        plt.axis('equal')

    def orbit3D(self):
        """Create 3D plot of orbits."""
        if not self.isgenerated:
            self.generateSystem()

        fig = plt.figure()  # creates figure
        ax = fig.add_subplot(111, projection='3d')  # create 3D subplot

        # running through celestial bodies:
        for i in range(self.numBods):
            if i == 0:
                plottype = "."
            else:
                plottype = "-"

            # plotting orbit of current body:
            ax.plot(self.bodyPos[:, 3*i],
                    self.bodyPos[:, 3*i + 1],
                    self.bodyPos[:, 3*i + 2],
                    plottype,
                    label=self.bodynames[i])

        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        ax.set_zlabel("z [AU]")
        ax.legend()

    def plotEnergy(self):
        if not self.isgenerated:
            self.generateSystem()

        totenergy = self.energy[:, 0] + self.energy[:, 1]

        plt.figure()
        plt.plot(self.times,
                 self.energy[:, 0],
                 label="Kinetic energy")
        plt.plot(self.times,
                 self.energy[:, 1],
                 label="Potential energy")
        plt.plot(self.times,
                 totenergy, '--',
                 label="Total energy")

        plt.xlabel("Time [year]")
        plt.ylabel("Energy")
        plt.legend()
        plt.grid()

    def plotAngMomMagnitude(self):
        if not self.isgenerated:
            self.generateSystem()

        angmommag = np.sqrt(self.angmom[:, 0]**2 +
                            self.angmom[:, 1]**2 +
                            self.angmom[:, 2]**2)

        plt.figure()
        plt.plot(self.times,
                 angmommag,
                 label="Angular momentum magnitude")

        plt.xlabel("Time [year]")
        plt.ylabel("Angular momentum")
        plt.legend()
        plt.grid()


bodynames = ["Sun", "Mercury", "Venus", "Earth", "Mars",
             "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

print("""Write se for Sun-Earth simulation,
sej for Sun-Earth-Jupiter,
sm for Sun-Mercury,
and ss for entire Solar System.""")

runflag = input("Choose run: ")
numTimesteps = int(eval(input("Number of time steps N = ")))
dt = float(eval(input("Size of time step dt = ")))

if runflag == "se":
    init_file = "earth-sun-init.txt"
    bodynames = [bodynames[0], bodynames[3]]

elif runflag == "sej":
    init_file = "sun-earth-jupiter-2020-Oct-19-00:00:00.init"
    bodynames = [bodynames[0], bodynames[3], bodynames[5]]

elif runflag == "sm":
    init_file = "sun-mercury-2020-Oct-19-00:00:00.init"
    bodynames = [bodynames[0], bodynames[1]]

elif runflag == "ss":
    init_file = "sun-and-friends-2020-Oct-19-00:00:00.init"


write_limit = 1
integration_method = "VelocityVerlet"
posfile = "positions.xyz"
momenfile = "energies.dat"

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

sun_earth.orbit3D()
sun_earth.orbit2D()
sun_earth.plotEnergy()
sun_earth.plotAngMomMagnitude()

plt.show()
# test_cpp()
clean()
