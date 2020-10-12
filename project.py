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


print("Set up run:")
numTimesteps = 1000
# dt = eval(input("Time step = "))
filename = "positions.xyz"

# Tstop = dt*numTimesteps

build_cpp()
run(["./main.exe", f"{numTimesteps}", filename], cwd=src)

with open(rootdir + "/data/" + filename) as infile:
    colsPrBod = 3
    header1 = infile.readline().split()
    numBods = int(header1[3][0])

    for i in range(2):
        infile.readline()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(numTimesteps):
        xyz = read_cords(infile, numBods, colsPrBod)
        if i > 0:
            for j in range(numBods):
                if j == 0:
                    plottype = "b."
                else:
                    plottype = "-"

                ax.plot([xyz[j, 0], old[j, 0]],
                        [xyz[j, 1], old[j, 1]],
                        [xyz[j, 2], old[j, 2]], plottype)

        old = xyz.copy()

plt.plot()
plt.show()


# test_cpp()
clean()
