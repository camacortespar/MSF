CXX = g++
CXXFLAGS = -std=c++17 -O3 -g
SPICE = -I/home/wind/naif/cspice/include /home/wind/naif/cspice/lib/cspice.a -lm
SANITFLAGS = -fsanitize=address -fsanitize=leak -fsanitize=undefined





all: %.x


%.x: %.cpp  
	$(CXX) $(CXXFLAGS) $< -o $@ $(SPICE)

clean:
	rm -rf *.x a.out *.x.*
