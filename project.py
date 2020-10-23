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


def init_file_maker(mat, filename):
    """
    Function writing init file from matrix containing initial conditions.
    Args:
        mat: numpy array of shape (n, 7) containing inital conditions
             for n bodies in the form: x0 y0 z0 vx0 vy0 vz0 m
        filename: string - name of file to create.
    """
    with open(rootdir + "/data/" + filename, "w") as outfile:
        for i in range(len(mat[:, 0])):
            for j in range(len(mat[0, :])):
                outfile.write(f"{mat[i, j]}")
                if j < len(mat[0, :])-1:
                    outfile.write(" ")
                else:
                    outfile.write("\n")


def test_cpp():
    """Function for running unit-tests."""
    run(["make", "test"], cwd=src)
    run("./test_main.exe", cwd=src)


def benchmark_cpp():
    """
    Function for running benchmark of euler and verlet integration using
    entire solar system.
    """
    # building and running benchmark program:
    run(["make", "benchmark"], cwd=src)
    run("./benchmark.exe", cwd=src)

    # reading benchmark data from file:
    benchmark_times = np.genfromtxt(rootdir +
                                    "/data/benchmarkdata.dat", skip_header=2)

    # taking mean and standard deviation:
    euler_mean = np.mean(benchmark_times[:, 0])
    euler_std = np.std(benchmark_times[:, 0])

    verlet_mean = np.mean(benchmark_times[:, 1])
    verlet_std = np.std(benchmark_times[:, 1])

    # printing results to terminal:
    header1 = "Time spent solving Solar system averaged over 50 runs."
    header2 = "N = 1000000, dt = 0.000248"
    eulerstr = f"Euler: {euler_mean:.4e} s \u00B1 {euler_std:.4e} s"
    verletstr = f"Verlet: {verlet_mean:.4e} s \u00B1 {verlet_std:.4e} s"
    print(header1)
    print(header2)
    print(eulerstr)
    print(verletstr)

    # writing results to file in data directory:
    with open(rootdir + "/data/benchmark_sun_and_friends.dat", "w") as output:
        output.write(header1 + "\n")
        output.write(header2 + "\n")
        output.write(eulerstr + "\n")
        output.write(verletstr + "\n")


class SolarSystem:
    """
    SolarSystemFiles is a class for reading data from SolarSystem C++
    program, and plotting orbits and energy/momentum.
    """

    def __init__(self, numTimesteps, dt, write_limit, integration_method,
                 init_file, posfile, momenfile, bodynames, correction, beta=2):
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
        self.write_limit = write_limit
        self.integration_method = integration_method
        # adding path to filenames:
        self.init_file = rootdir + "/data/" + init_file
        self.posfile = rootdir + "/data/" + f"{beta*10:.0f}" + posfile
        self.momenfile = rootdir + "/data/" + f"{beta*10:.0f}" + momenfile

        self.bodynames = bodynames
        self.correction = correction
        self.beta = beta
        self.runflag = "normal"

        # ensuring that maximum 1000 lines of data are read from each file,
        # this is to limit memory use.
        if (self.numTimesteps//write_limit > 10000):
            self.everyNlines = self.numTimesteps // (write_limit*10000)
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
                f"{self.numTimesteps}",
                f"{self.dt}",
                f"{self.write_limit}",
                self.integration_method,
                self.init_file,
                self.posfile,
                self.momenfile,
                self.correction,
                f"{self.beta}",
                self.runflag
            ],
            cwd=src
        )

        self.readData()  # reading data from files
        print("\nDone generating data")
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

        if self.runflag != "sm":
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
                    # self.energy is an Nx2-array containing the Kinetic energy
                    # each time step as the first column and potential energy
                    # each time step as the second column.

    def orbit2D(self, number_of_bodies=None, center_on_sun=False,
                axis='equal'):
        """
        Method for creating 2D plot of orbits,
        by plotting the x and y coordinates of the bodies' positions.
        Args:
            number_of_bodies: int - number of bodies to plot,
                                    counting from zero.
        """
        if not self.isgenerated:  # generating data if necessary.
            self.generateSystem()
        if number_of_bodies is None:
            number_of_bodies = self.numBods
        if center_on_sun:
            self.moveToSunFrame()

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
            f", Simulated time = {self.dt*self.numTimesteps:.2f} years\n"
        )
        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        plt.legend()
        plt.grid()
        plt.axis(axis)

    def orbit3D(self, number_of_bodies=None, center_on_sun=False):
        """
        Create 3D plot of orbits.
        Args:
            number_of_bodies: int - number of bodies to plot,
                                    counting from zero.
        """
        if not self.isgenerated:  # generating data if necessary.
            self.generateSystem()
        if number_of_bodies is None:
            number_of_bodies = self.numBods
        if center_on_sun:
            self.moveToSunFrame()

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
        """
        Plotting kinetic, potential and total energy for system as
        function of time steps. Values are scaled by maximum magnitude of
        total energy.
        """
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
            f", Simulated time = {self.dt*self.numTimesteps:.2f} years\n"
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
            f", Simulated time = {self.dt*self.numTimesteps:.2f} years\n"
        )
        plt.xlabel("N [number of time steps]")
        plt.ylabel("Angular momentum $[L_{tot}/max(|L_{tot}|)]$")
        plt.legend()
        plt.grid()

    def moveToSunFrame(self):
        """Method for moving positions to Sun's frame of refrence."""
        if not self.isgenerated:
            self.generateSystem()

        for i in range(self.numBods):
            self.bodyPos[:, 3*i:3*i+3] = (self.bodyPos[:, 3*i:3*i+3] -
                                          self.bodyPos[:, :3])
        # self.bodyPos[:, :3] = np.zeros(
        #     (self.numTimesteps//(self.everyNlines*self.write_limit), 3)
        # )

    def perihelionAngle(self):
        """Method for calculating all perihelionAngles in array of body 1."""
        self.everyNlines = 1
        self.write_limit = 1
        self.runflag = "sm"
        if not self.isgenerated:
            self.generateSystem()
        if self.correction == "nonrel":
            self.moveToSunFrame()

        # finding perihelion:
        # T = 0.2411
        # N = (self.numTimesteps - int(T/self.dt))//(self.everyNlines *
        #                                            self.write_limit)
        # print(N-1)
        # print(self.bodyPos.shape)
        # rvec = self.bodyPos[N:, 3:]
        # print(rvec.shape)
        # r = np.sqrt(rvec[:, 0]**2 + rvec[:, 1]**2 + rvec[:, 2]**2)
        # self.a = np.argmin(r)

        self.xp = self.bodyPos[:, 3]
        self.yp = self.bodyPos[:, 4]
        self.thetaP = np.arctan2(self.yp, self.xp)*206264.806


# name of bodies used in project:
bodynames = ["Sun", "Mercury", "Venus", "Earth", "Mars",
             "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
beta = 2
runflag = "start"
betaflag = "n"
write_limit = 1
center_on_sun = True
axis = 'equal'
num_bods_to_plot = None


# asking for imput:
print("""Write se for Sun-Earth simulation,
sej for Sun-Earth-Jupiter,
sm for Sun-Mercury,
and ss for entire Solar System.
Write test to run unit-tests,
or b for benchmark.""")


while (runflag != "se" and runflag != "sej" and runflag != "sm" and
       runflag != "ss" and runflag != "test" and runflag != "b"):

    runflag = input("Choose run: ")


if runflag != "test" and runflag != "b":
    # asking for imput parameters from user:
    numTimesteps = int(eval(input("\nNumber of time steps N = ")))
    dt = float(eval(input("Size of time step dt = ")))
    if runflag != "sm":
        limit_write = input(
            "\nOnly write the data from 10000 evenly spaced time steps? y/n: "
        )
        if limit_write == "y" and numTimesteps > 10000:
            write_limit = numTimesteps//10000

if runflag == "se":  # initial data for sun_earth run:
    init_file = "sun-earth.init"
    bodynames = [bodynames[0], bodynames[3]]

    escvelflag = input("Run escape velocity test? y/n: ")
    if escvelflag == "y":
        escvel = float(eval(input("Enter initial escape velocity: ")))

        # creating file with initial conditions:
        init_file = "sun-earth-escvel.init"
        mat = np.zeros((2, 7))
        mat[0, 6] = 1  # sun mass
        mat[1, 0] = 1  # earth x0
        mat[1, 4] = escvel  # earth vx0
        mat[1, 6] = 3e-6  # earth mass
        init_file_maker(mat, init_file)

        axis = "auto"

    else:
        betaflag = input("Run sun earht with varying beta? y/n: ")
        if betaflag == "y":
            beta_array = np.linspace(2, 3, 4)

            initflag = input("Use eliptical orbit? y/n: ")
            if initflag == "y":
                init_file = "sun-earth-peturbed.init"

elif runflag == "sej":  # initial data for sun_earth_jupiter run:
    init_file = "sun-earth-jupiter-2020-Oct-19-00:00:00.init"
    bodynames = [bodynames[0], bodynames[3], bodynames[5]]
    jup_m_fac = [1, 10, 1000]
    init_files = [init_file]

    for i in jup_m_fac:
        if i > 1:
            initconditions = np.genfromtxt(rootdir + "/data/" + init_file)
            initconditions[2, 6] *= i
            file = f"sej_{i}" + init_file
            init_file_maker(initconditions, file)
            init_files.append(file)


elif runflag == "sm":  # initial data for sun_mercury run:
    init_file = "sun-mercury.init"
    bodynames = [bodynames[0], bodynames[1]]


elif runflag == "ss":  # initial data for entire SolarSystem run:
    center_on_sun = False
    init_file = "sun-and-friends-2020-Oct-19-00:00:00.init"

    innerboddies = input("Only plot sun + 5 innermost planets? y/n: ")
    if innerboddies == "y":
        num_bods_to_plot = 6

if (runflag != "test") and (runflag != "b"):  # setting up run:

    integration_method = input("""
Choose integration method:
    Write fe for forward Euler,
    or vv for Velocity-Verlet:
    """)

    # setting up algorithm and output filenames:
    if integration_method == "vv":
        integration_method = "VelocityVerlet"
    elif integration_method == "fe":
        integration_method = "Euler"

    posfile = runflag + "_" + integration_method + "_" + "positions.xyz"
    momenfile = runflag + "_" + integration_method + "_" + "energies.dat"
    # initalising instance of SolarSystem class with parameters:
    if runflag != "sm" and betaflag == "y":
        posdict = {}
        for beta in beta_array:
            system = SolarSystem(
                numTimesteps,
                dt,
                write_limit,
                integration_method,
                init_file,
                posfile,
                momenfile,
                bodynames,
                "nonrel",
                beta
            )

            system.orbit3D(center_on_sun=True)
            system.orbit2D(axis=axis)
            system.plotEnergy()
            system.plotAngMomMagnitude()
            plt.show()

    elif runflag != "sm" and runflag != "sej" and betaflag != "y":
        system = SolarSystem(
            numTimesteps,
            dt,
            write_limit,
            integration_method,
            init_file,
            posfile,
            momenfile,
            bodynames,
            "nonrel",
            beta
        )

    elif runflag == "sej":
        for i in range(len(init_files)):
            system = SolarSystem(
                numTimesteps,
                dt,
                write_limit,
                integration_method,
                init_files[i],
                posfile,
                momenfile,
                bodynames,
                "nonrel",
                beta
            )
            print(f"\nPlotting for Jupiter mass multiplier = {jup_m_fac[i]}x")
            system.orbit3D(number_of_bodies=num_bods_to_plot,
                           center_on_sun=center_on_sun)
            system.orbit2D(number_of_bodies=num_bods_to_plot,
                           axis=axis)
            system.plotEnergy()
            system.plotAngMomMagnitude()

            plt.show()

    if runflag != "sm" and runflag != "sej" and betaflag != "y":
        # calling the various plot functions:
        system.orbit3D(number_of_bodies=num_bods_to_plot,
                       center_on_sun=center_on_sun)
        system.orbit2D(number_of_bodies=num_bods_to_plot,
                       axis=axis)
        system.plotEnergy()
        system.plotAngMomMagnitude()

        plt.show()

    elif runflag == "sm":
        system = SolarSystem(
            numTimesteps,
            dt,
            write_limit,
            integration_method,
            init_file,
            posfile,
            momenfile,
            bodynames,
            "nonrel",
            beta
        )

        system2 = SolarSystem(numTimesteps, dt, write_limit,
                              integration_method, init_file,
                              "rel_" + posfile, "rel_" + momenfile,
                              bodynames, "rel")
        system.perihelionAngle()
        system2.perihelionAngle()
        # system.orbit2D()
        # system2.orbit2D()

        plt.show()

        thetaPnonrel = system.thetaP[-1]
        thetaPrel = system2.thetaP[-1]

        print(f"ThetaP without correction = {thetaPnonrel:g}''")
        print(f"ThetaP with  correction = {thetaPrel:g}''")
        # print(f"{system.thetaP[-1]:e}", f"{system2.thetaP[-1]:e}")


elif runflag == "test":
    test_cpp()

elif runflag == "b":
    benchmark_cpp()
