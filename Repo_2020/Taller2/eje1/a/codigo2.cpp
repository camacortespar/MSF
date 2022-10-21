#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>



struct Particle{

  int position[2] = {0,0};

  void Move(int step, int direction, const CONFIG &config);
  int Getcell(const CONFIG &config);
};

int Particle::Getcell(const CONFIG &config){

    int X = (position[0] + config.latticesize/2)*config.gridsize/config.latticesize;
    int Y = (position[1] + config.latticesize/2)*config.gridsize/config.latticesize;

    return X + Y*config.gridsize;
}

void Particle::Move(int step, int direction, const CONFIG &config){

  if (labs(position[direction] + step) != config.latticesize/2 + (1 - step)/2){

        position[direction] += step;

    }

}
void start( Vec_i &Cells, Vec_p &Particles);

void time_step( int random_particle, int step, int direction, Vec_p &Particles);

double dropsize(const Vec_p &Particles);


int main(void)
{

  CONFIG config;
  config.read("Data/init_data.txt");

  Vec_p Particles(config.nmolecules);
  Vec_i Cells(config.gridsize*config.gridsize,0);

  start(config, Cells, Particles);

  std::mt19937 gen(config.seed);
  std::uniform_int_distribution<int> dis_move(0, 1);
  std::uniform_int_distribution<int> dis_particle(0,config.nmolecules-1);

  double Size=0;
  int random_particle = 0, step = 0, direction = 0;

  for(int t = 0; t <= config.tmax; t++ ){

    random_particle = dis_particle(gen);   //escoge una particula al azar                                                                                           
    step = dis_move(gen)*2 - 1;           //genera un numero aleatorio 1 o -1 (1: arriba o derecha -1:abajo o izquierda)                                         
    direction = dis_move(gen);           //genera un numero aleatorio 0 o 1 (0 para x 1 para y)                                                  

    time_step(config, random_particle, step, direction, Cells, Particles);

    if (t%config.resolution == 0){

      Size= dropsize(config,Particles);

      std::cout << t <<"\t" << Size << std::endl;
    }
  }
  return 0;
}

void time_step( int random_particle, int step, int direction, Vec_p &Particles){

    int m = 0;

    m = Particles[random_particle].Getcell(config);
    Cells[m] -= 1;
    Particles[random_particle].Move(step,direction,config);
    m = Particles[random_particle].Getcell(config);
    Cells[m] += 1;
}







void start( Vec_i &Cells, Vec_p &Particles)
{
    int n = std::sqrt(nmolecules);
    int m = 0;

    for(int i = 0; i < nmolecules; i++){
        //inicializa las particulas en un cuadrado                              
        Particles[i].position[0] = i%n - n/2;
        Particles[i].position[1] = i/n - n/2;
        //calcula las particulas en cada celda                                  
        m = Particles[i].Getcell(config);
        Cells[m] += 1;
    }

}

double dropsize(const CONFIG &config,const Vec_p &Particles)
{
  double suma=0,x=0,y=0,r=0,size=0;

    for(auto i: Particles){
      x=i.position[0]-0.5;
      y=i.position[1]-0.5;
      r=std::hypot(x,y);
      suma +=std::pow(r,2);
  }
    size=std::sqrt(suma/config.nmolecules);

    return size;
}

void time_step( int random_particle, int step, int direction, Vec_i &Cells, Vec_p &Particles){

    int m = 0;

    Particles[random_particle].Move(step,direction,config);
}
double dropsize(const Vec_p &Particles)
{
  double suma=0,x=0,y=0,r=0,size=0;
  
    for(auto i: Particles){
      x=i.position[0]-0.5;
      y=i.position[1]-0.5;
      r=std::hypot(x,y);
      suma +=std::pow(r,2);
  }
    size=suma/nmolecules;

    return size;
}
