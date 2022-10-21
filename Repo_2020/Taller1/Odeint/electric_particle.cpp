//

#include <iostream>
#include <array>
#include <cmath>

#include <boost/numeric/odeint.hpp>

//using namespace boost::numeric::odeint;

const double Bx = 0;
const double By = 0;
const double Bz = 0;
const double Ex = 0;
const double Ey = 0;
const double Ez =-5e30;
const double K=4.65e57;
const double e=-1;
const double tmin=0;
const double tmax=10;
const double dt=0.1;


typedef std::array< double , 6 > state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
  dxdt[0] = x[3];
  dxdt[1] = x[4];
  dxdt[2] = x[5];
  
  double r=std::hypot(x[1],x[2],x[3]);
  double r3=std::pow(r,3);
    if(r3==0){r3=0.0001;}
  
  dxdt[3] = e*(K*x[0]/r3+Ex+x[4]*Bz-x[5]*By);
  dxdt[4] = e*(K*x[1]/r3+Ey+x[5]*Bx-x[3]*Bz);
  dxdt[5] = e*(K*x[2]/r3+Ez+x[3]*By-x[4]*Bx);
}

void write_lorenz( const state_type &x , const double t )
{
  double rho=std::hypot(x[0],x[1]);
  std::cout << rho << '\t' << x[2] << std::endl;
}

int main(int argc, char **argv)
{
  state_type x = { 1,0,0,10,0,0}; // initial conditions
  
                                     
    boost::numeric::odeint::integrate( lorenz , x , tmin , tmax , dt , write_lorenz );
}
