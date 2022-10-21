#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include "Random64.h"

const int nmolecules = 2400;
const  int gridsize = 256;
const int latticesize = 256;
const int tmax = 350;
const  int seed = 0;
const  int resolution = 5;
const double p=0.25;
const double p0=0.5;

class Particle{
private:
public:
  int position[2] = {0,0};  

  void Move(int step, int direction);
};

void Particle::Move(int step, int direction){

  if (labs(position[direction] + step) != latticesize/2 + (1 - step)/2){

        position[direction] += step;
    }
}
//n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i];

typedef std::vector<Particle> Vec_p;
typedef std::vector<int> Vec_i;
typedef std::vector<double> Vec_d;

void start( Vec_p &Particles);

double dropsize(const Vec_p &Particles);

int main(void)
{

  Vec_p Particles(nmolecules);
  
  start(Particles);

  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dis_move(0, 1);
  
  double Size=0;
  int step = 0, direction = 0;
  double probabilidad=0;
  
  for(int t = 0; t <= tmax; t++ ){
    for(int i=0; i<nmolecules; i++){
      
    probabilidad=dis_move(gen);
    
    if(probabilidad<p0)
      {step=1;direction=0;}
    else if(p0<=probabilidad<p)
      {step=1;direction=1;}
    else if(p<=probabilidad<1-2*p-p0)
      {step=-1;direction=0;}
    else if(1-2*p-p0<=probabilidad)
      {step=-1;direction=1;}
 
    Particles[i].Move(step,direction);
    }
    if (t%resolution == 0){
        
      Size= dropsize(Particles);
      
      std::cout << t <<"\t" << Size << std::endl;
    }
  }  
  return 0;
}

void start( Vec_p &Particles)
{
  Crandom ran64(2);
  int n = std::sqrt(nmolecules);
  int m = 0;
  int Lx=latticesize;
    
  for(int i = 0; i < nmolecules; i++){
    //inicializa las particulas
    int ix=(int) ran64.gauss(0,16); if(ix<0) ix=0; if(ix>(Lx-1)) ix=Lx-1;
    int iy=(int) ran64.gauss(0,16); if(ix<0) ix=0; if(ix>(Lx-1)) iy=Lx-1;
    Particles[i].position[0] = ix;
    Particles[i].position[1] = iy;
    }
    
}


double dropsize(const Vec_p &Particles)
{
  double suma=0,x=0,y=0,r=0,size=0,var=0;

    for(auto i: Particles){
      x=i.position[0]-16;
      y=i.position[1]-16;
      r=std::hypot(x,y);
      suma +=r;
  }
    size=suma/nmolecules;

    for(auto i: Particles)
      var=+std::pow(std::hypot(i.position[0]-16,i.position[1]-16)-size,2);
      
    return size*size;
}
