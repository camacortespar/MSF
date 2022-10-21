CXX = g++
CXXFLAGS = -std=c++17 -O3 -g
SANITFLAGS = -fsanitize=address -fsanitize=leak -fsanitize=undefined





all: main.x


main.x: CA_Difusion_1D_continuo.cpp 
	$(CXX) $(CXXFLAGS) $(SANITFLAGS) $<  


clean:
	rm -rf *.x a.out *.x.*
