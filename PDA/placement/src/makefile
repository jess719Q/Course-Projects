all: main.cpp detailedPM.cpp detailedPM.h
	g++ -std=gnu++17 -O3 -fopenmp -march=native -funroll-loops -DNDEBUG -o ../bin/hw3 main.cpp detailedPM.cpp
clean:
	rm -f ../bin/hw3