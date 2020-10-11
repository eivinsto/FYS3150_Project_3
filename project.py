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
pwd = os.getcwd()
wd = pwd + "/src"


def build_cpp():
    """Function for building c++ program."""
    run(["make", "all"], cwd=wd)


def test_cpp():
    """Function for running unit-tests."""
    run(["make", "test"], cwd=wd)
    run("./test_main.exe", cwd=wd)


def clean():
    """Function for cleaning datafiles in src directory."""
    run(["make", "cleandat"], cwd=wd)
