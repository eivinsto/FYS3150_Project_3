CC := g++
FLAGS := -Wall -O2 -march=native
LIBS := -larmadillo
OFILES := solar_integrator.o solar_system.o celestial_body.o

.PHONY: all clean cleanbin test

# compile code:
all: main.exe

# compile program:
main.exe: main.o $(OFILES)
	$(CC) $(FLAGS) main.o $(OFILES) -o main.exe $(LIBS)

main.o: main.cpp
	$(CC) $(FLAGS) main.cpp -c

solar_integrator.o: solar_integrator.cpp
	$(CC) $(FLAGS) solar_integrator.cpp -c

solar_system.o: solar_system.cpp
	$(CC) $(FLAGS) solar_system.cpp -c

celestial_body.o: celestial_body.cpp
	$(CC) $(FLAGS) celestial_body.cpp -c

# project cleaning:
clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.dat

# data cleaning:
cleandat:
	rm -f *.dat

print-% : ; @echo $* = $($*)
