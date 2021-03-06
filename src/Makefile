CC := g++
FLAGS := -Wall -Wextra -O2 -march=native
DEBUG := -g
LIBS := -larmadillo
OFILES := solar_integrator.o solar_system.o celestial_body.o

.PHONY: all clean benchmark test

# compile code:
all: main.exe

# compile program:
main.exe: main.o $(OFILES)
	$(CC) $(FLAGS) $(DEBUG) main.o $(OFILES) -o main.exe $(LIBS)

main.o: main.cpp
	$(CC) $(FLAGS) $(DEBUG) main.cpp -c

solar_integrator.o: solar_integrator.cpp
	$(CC) $(FLAGS) $(DEBUG) solar_integrator.cpp -c

solar_system.o: solar_system.cpp
	$(CC) $(FLAGS) $(DEBUG) solar_system.cpp -c

celestial_body.o: celestial_body.cpp
	$(CC) $(FLAGS) $(DEBUG) celestial_body.cpp -c


# unittest program:
test: test_main.exe

test_main.exe: test_main.o test_functions.o $(OFILES)
	$(CC) $(FLAGS) test_main.o test_functions.o $(OFILES) -o test_main.exe $(LIBS)

test_main.o: test_main.cpp
	$(CC) $(FLAGS) test_main.cpp -c

test_functions.o: test_functions.cpp
	$(CC) $(FLAGS) test_functions.cpp -c

# benchmark
benchmark: benchmark.exe

benchmark.exe: benchmark.o $(OFILES)
	$(CC) $(FLAGS) benchmark.o $(OFILES) -o benchmark.exe $(LIBS)

benchmark.o: benchmark.cpp
	$(CC) $(FLAGS) benchmark.cpp -c

# project cleaning:
clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.dat

print-% : ; @echo $* = $($*)
