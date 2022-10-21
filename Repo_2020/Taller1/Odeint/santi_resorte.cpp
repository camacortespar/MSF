//

#include <iostream>
#include <array>
#include <cmath>

#include <boost/numeric/odeint.hpp>

//using namespace boost::numeric::odeint;

const double k = 1;
const double m = 10;
const double l=1;
const double d=0.2;
const double tmin=0;
const double tmax=30;
const double dt=0.00001;


typedef std::array< double , 4 > state_type;

void lorenz( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = x[1];
    dxdt[1] = -k*(x[0]+l)/m-d*x[1];
    dxdt[2]=x[3];
    dxdt[3] =2*x[1]*x[3]/(x[0]+l);
}
void write_lorenz( const state_type &x , const double t )
{
  double a=(x[0]+l)*std::cos(x[2]);
  double b=(x[0]+l)*std::sin(x[2]);
  
  std::cout << '\t' << a << '\t' << b << std::endl;
}

int main(int argc, char **argv)
{
  state_type x = { 1.5 , 0 , 0, 0.1 }; // initial conditions
                                     
  boost::numeric::odeint::integrate( lorenz , x , tmin , tmax , dt , write_lorenz );
}
