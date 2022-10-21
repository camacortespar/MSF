//

#include <iostream>
#include <array>
#include <cmath>

#include <boost/numeric/odeint.hpp>

//using namespace boost::numeric::odeint;

const double F = 0;
const double P = 200;
const double tmin=0;
const double tmax=100;
const double dt=0.1;


typedef std::array< double , 4 > state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = x[2];
    dxdt[1] = x[3];
    if(x[0]==1){
      dxdt[2]=1;}
    else{
      dxdt[2] = -P*P/std::pow(x[0],3)+x[0]/(std::pow(std::hypot(x[0],x[1]),3));}
    if(x[0]==0 and x[1]==0){
      dxdt[3]=0;}
    else{
      dxdt[3] =F+x[1]/(std::pow(std::hypot(x[0],x[1]),3));}
}

void write_lorenz( const state_type &x , const double t )
{
  std::cout << t << '\t' << x[0] << '\t' << x[1] << std::endl;
}

int main(int argc, char **argv)
{
  state_type x = { 1 , -2 , -1, -1 }; // initial conditions
  
                                     
    boost::numeric::odeint::integrate( lorenz , x , tmin , tmax , dt , write_lorenz );
}
