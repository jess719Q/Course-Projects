all: main.cpp router.cpp router.h
	g++ -std=gnu++17 -O3 -fopenmp -march=native -funroll-loops -DNDEBUG -o ../bin/hw4 main.cpp router.cpp
# 	g++ -std=gnu++17 -o ../bin/hw4 main.cpp router.cpp
clean:
	rm -f ../bin/hw4