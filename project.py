"""
Python script to interface with project code.
"""
from subprocess import run
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import shutil

# retriveing working directories:
rootdir = os.getcwd()
src = rootdir + "/src"


def main():
    build_cpp()
    test_cpp()
    clean()


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


if __name__ == '__main__':
    main()
