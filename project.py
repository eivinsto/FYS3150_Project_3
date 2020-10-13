"""
Python script to interface with project code.
"""
from subprocess import run
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import shutil

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


def read_cords(infile, numBods, colsPrBod):
    xyz = np.zeros((numBods, 3))
    line = infile.readline().split()
    for i in range(numBods):
        for j in range(colsPrBod):
            xyz[i, j] = float(line[colsPrBod*i + j])
    return xyz


def clean(files="dat"):
    """Function for cleaning datafiles in src directory."""
    if files == "dat":
        run(["make", "cleandat"], cwd=src)
    if files == "all":
        run(["make", "clean"], cwd=src)


class SolarSystemFiles:
    """SolarSystemFiles is a class for reading data from SolarSystem C++
    program, and plotting orbits and energy/momentum."""

    def __init__(self, filename, bodynames):
        """Class constructor.
        Args:
            filename: string - Name of datafile containing plannet positions.
            bodynames: list of strings - List of names of each celestial body.
        """
        self.filename = rootdir + "/data/" + filename
        self.bodynames = bodynames

    def readHeader(self):
        """Read header from datafile and set number of bodies."""
        with open(self.filename) as infile:
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
            self.filename,
            skip_header=3,
            usecols=np.arange(
                self.colsPrBod*currentbod,
                self.colsPrBod*currentbod + self.colsPrBod
            )
        )

    def orbit3D(self):
        """Create 3D plot of orbits."""
        fig = plt.figure()  # creates figure
        ax = fig.add_subplot(111, projection='3d')  # create 3D subplot
        self.readHeader()  # reading data from header

        # running through celestial bodies:
        for i in range(self.numBods):
            if i == 0:
                plottype = "."
            else:
                plottype = "-"

            self.readBodyData(i)  # reading position data of current body
            # plotting orbit of current body:
            ax.plot(self.bodyPos[:, 0],
                    self.bodyPos[:, 1],
                    self.bodyPos[:, 2],
                    plottype,
                    label=self.bodynames[i])

        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        ax.set_zlabel("z [AU]")
        ax.legend()
        plt.show()


print("Set up run:")
numTimesteps = 1000
# dt = eval(input("Time step = "))
filename = "positions.xyz"
bodynames = ["Sun", "Earth"]
sun_earth = SolarSystemFiles(filename, bodynames)
# Tstop = dt*numTimesteps

build_cpp()
run(["./main.exe", f"{numTimesteps}", filename], cwd=src)
sun_earth.orbit3D()


# test_cpp()
clean()
