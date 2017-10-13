2pcfv1 : tools.o 2pcfv1.o
	g++ 2pcfv1.o tools.o -o 2pcfv1 -std=c++11

tools.o : tools.h tools.cpp
	g++ -c tools.cpp -o tools.o -std=c++11

2pcfv1.o : tools.h 2pcfv1.cpp
	g++ -c 2pcfv1.cpp -o 2pcfv1.o
