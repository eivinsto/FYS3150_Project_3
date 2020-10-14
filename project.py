"""
Python script to interface with project code.
"""
from subprocess import run
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import sys

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


class SolarSystemFiles:
    """SolarSystemFiles is a class for reading data from SolarSystem C++
    program, and plotting orbits and energy/momentum."""

    def __init__(self, posfile, momenfile, bodynames, numTimesteps, dt):
        """Class constructor.
        Args:
            posfile: string - Name of datafile containing plannet positions.
            momenfile: string - Name of datafile containing energy and angular
                                momentum.
            bodynames: list of strings - List of names of each celestial body.
            numTimesteps: integer - number of time steps to take.
            dt: float (64bit) - size of time step.
        """
        self.posfile = rootdir + "/data/" + posfile
        self.momenfile = rootdir + "/data/" + momenfile
        self.bodynames = bodynames
        self.numTimesteps = numTimesteps
        self.dt = dt

    def readHeader(self):
        """Read header from datafile and set number of bodies."""
        with open(self.posfile) as infile:
            self.colsPrBod = 3
            header1 = infile.readline().split()
            self.numBods = int(header1[3][0])

    def readBodyData(self, currentbod):
        """Reads position-data for current body from datafile.
        Result is stored as x y z coulmns in 3xN array self.bodyPos.
        Arg:
            currentbod: Integer - index of current celestial body.
        """
        self.bodyPos = np.genfromtxt(
            self.posfile,
            skip_header=3,
            usecols=np.arange(
                self.colsPrBod*currentbod,
                self.colsPrBod*currentbod + self.colsPrBod
            )
        )

    def orbit2D(self):
        """Create 2D plot of orbits."""
        self.readHeader()  # reading data from header
        plt.figure()  # creates figure
        # limiting number of datapoints plotted
        if self.numTimesteps > 1000:
            plotStep = self.numTimesteps//1000
        else:
            plotStep = 1

        # running through celestial bodies:
        for i in range(self.numBods):
            if i == 0:
                plottype = "."
            else:
                plottype = "-"
            self.readBodyData(i)
            plt.plot(self.bodyPos[::plotStep, 0],
                     self.bodyPos[::plotStep, 1],
                     plottype,
                     label=self.bodynames[i])

        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        plt.legend()
        plt.grid()
        plt.axis('equal')

    def orbit3D(self):
        """Create 3D plot of orbits."""
        fig = plt.figure()  # creates figure
        ax = fig.add_subplot(111, projection='3d')  # create 3D subplot
        self.readHeader()  # reading data from header
        # limiting number of datapoints plotted:
        if self.numTimesteps > 1000:
            plotStep = self.numTimesteps//1000
        else:
            plotStep = 1

        # running through celestial bodies:
        for i in range(self.numBods):
            if i == 0:
                plottype = "."
            else:
                plottype = "-"

            self.readBodyData(i)  # reading position data of current body
            # plotting orbit of current body:
            ax.plot(self.bodyPos[::plotStep, 0],
                    self.bodyPos[::plotStep, 1],
                    self.bodyPos[::plotStep, 2],
                    plottype,
                    label=self.bodynames[i])

        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        ax.set_zlabel("z [AU]")
        ax.legend()

    def plotEnergy(self):
        self.readHeader()
        if self.numTimesteps > 1000:
            plotStep = self.numTimesteps//1000
        else:
            plotStep = 1

        energy = np.genfromtxt(self.momenfile, skip_header=1, usecols=[0, 1])
        totenergy = energy[:, 0] + energy[:, 0]
        times = np.linspace(0, self.numTimesteps*dt, self.numTimesteps)

        plt.figure()

        plt.plot(times[::plotStep],
                 energy[::plotStep, 0],
                 label="Kinetic energy")
        plt.plot(times[::plotStep],
                 energy[::plotStep, 1],
                 label="Potential energy")
        plt.plot(times[::plotStep],
                 totenergy[::plotStep], '--',
                 label="Total energy")

        plt.xlabel("Time [year]")
        plt.ylabel("Energy")
        plt.legend()
        plt.grid()

    def plotAngMomMagnitude(self):
        self.readHeader()
        if self.numTimesteps > 1000:
            plotStep = self.numTimesteps//1000
        else:
            plotStep = 1

        angmom = np.genfromtxt(self.momenfile,
                               skip_header=1,
                               usecols=[2, 3, 4])
        angmommag = np.sqrt(angmom[:, 0]**2 +
                            angmom[:, 1]**2 +
                            angmom[:, 2]**2)
        times = np.linspace(0, self.numTimesteps*dt, self.numTimesteps)

        plt.figure()
        plt.plot(times[::plotStep],
                 angmommag[::plotStep],
                 label="Angular momentum magnitude")

        plt.xlabel("Time [year]")
        plt.ylabel("Angular momentum")
        plt.legend()
        plt.grid()


print("Set up run:")
numTimesteps = 1000
if len(sys.argv) >= 2:
    numTimesteps = int(eval(sys.argv[1]))
if len(sys.argv) >= 3:
    dt = float(eval(sys.argv[2]))
# dt = eval(input("Time step = "))
posfile = "positions.xyz"
momenfile = "energies.dat"
bodynames = ["Sun", "Earth"]
sun_earth = SolarSystemFiles(posfile, momenfile, bodynames, numTimesteps, dt)
# Tstop = dt*numTimesteps

build_cpp()
run(["./main.exe", f"{numTimesteps}", f"{dt}"], cwd=src)
print("Integration done!")
sun_earth.orbit3D()
sun_earth.orbit2D()
sun_earth.plotEnergy()
sun_earth.plotAngMomMagnitude()

plt.show()
# test_cpp()
clean()
