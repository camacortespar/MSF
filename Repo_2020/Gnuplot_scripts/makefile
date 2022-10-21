

all: main.x

main.x: atomo_hidrogeno.cpp
	g++ -O3 atomo_hidrogeno.cpp -o $@
	./main.x
	gnuplot plot_atomo-espacio.gp
	gnuplot plot_atomo-plano.gp


