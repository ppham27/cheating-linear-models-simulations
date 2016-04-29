CPPFLAGS=-stdlib=libc++ -std=c++0x

default: independentSimulation powerSimulation calculateVariance

simulate: simulate.o
	g++ $(CPPFLAGS) -c simulate.cpp

independentSimulation: simulate
	g++ $(CPPFLAGS) -o independentSimulation independentSimulation.cpp simulate.o -larmadillo

powerSimulation: simulate
	g++ $(CPPFLAGS) -o powerSimulation powerSimulation.cpp simulate.o -larmadillo

calculateVariance: simulate
	g++ $(CPPFLAGS) -o calculateVariance calculateVariance.cpp simulate.o -larmadillo

testBalance: simulate
	g++ $(CPPFLAGS) -o testBalance testBalance.cpp simulate.o -larmadillo
