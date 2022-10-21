CXX_P = g++
CXXFLAGS = -std=c++17 -O3 -g
DEPENDENCIES = eje1.cpp


all: eje1 eje2 eje3 eje4 animacion1 animacion3

eje1: eje1.cpp
	$(CXX_P) $(CXXFLAGS) eje1.cpp -o eje1.x
	./eje1.x
	gnuplot eje1.gp

eje2: eje2.cpp
	$(CXX_P) $(CXXFLAGS) eje2.cpp -o eje2.x
	./eje2.x
	gnuplot eje2.gp

eje3: eje3.cpp
	$(CXX_P) $(CXXFLAGS) eje3.cpp -o eje3.x
	./eje3.x
	gnuplot eje3.gp

eje4: eje4.cpp
	$(CXX_P) $(CXXFLAGS) eje4.cpp -o eje4.x
	./eje4.x
	gnuplot eje4.gp

animacion1: eje1_animacion.cpp
	$(CXX_P) $(CXXFLAGS) eje1_animacion.cpp -o eje1a.x
	./eje1a.x | gnuplot

animacion3: eje3_animacion.cpp
	$(CXX_P) $(CXXFLAGS) eje3_animacion.cpp -o eje3a.x
	./eje3a.x | gnuplot


clean:
	rm -rf *.x a.out *.x.* *.dat

