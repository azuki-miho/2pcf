2pcfv1 : tools.o 2pcfv1.o
	g++ 2pcfv1.o tools.o -o 2pcfv1 -std=c++11

2pcfv2 : tools.o 2pcfv2.o
	g++ 2pcfv2.o tools.o -o 2pcfv2 -std=c++11

2pcfv2_1 : tools.o 2pcfv2_1.o
	g++ -O3 2pcfv2_1.o tools.o -o 2pcfv2_1 -std=c++11

tools.o : tools.h tools.cpp
	g++ -O3 -c tools.cpp -o tools.o -std=c++11

2pcfv1.o : tools.h 2pcfv1.cpp
	g++ -c 2pcfv1.cpp -o 2pcfv1.o -std=c++11

2pcfv2.o : tools.h 2pcfv2.cpp
	g++ -c 2pcfv2.cpp -o 2pcfv2.o -std=c++11

2pcfv2_1.o : tools.h 2pcfv2_1.cpp
	g++ -O3 -c 2pcfv2_1.cpp -o 2pcfv2_1.o -std=c++11

selectgalaxy : tools.o selectgalaxy.o
	g++ selectgalaxy.o tools.o -o selectgalaxy -std=c++11

selectgalaxy.o : tools.h selectgalaxy.cpp
	g++ -c selectgalaxy.cpp -o selectgalaxy.o

generaterandom : tools.o generaterandom.o
	g++ generaterandom.o tools.o -o generaterandom -std=c++11

generaterandom.o : tools.h generaterandom.cpp
	g++ -c generaterandom.cpp -o generaterandom.o

log_linear : log_linear.o tools.o
	g++ log_linear.o tools.o -o log_linear -std=c++11

log_linear.o : tools.h log_linear.cpp
	g++ -c log_linear.cpp -o log_linear.o

ztor : tools.o ztor.o
	g++ ztor.o tools.o -o ztor -std=c++11

ztor.o : tools.h ztor.cpp
	g++ -c ztor.cpp -o ztor.o
