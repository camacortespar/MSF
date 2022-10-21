#include <iostream>
#include <array>
#include <boost/numeric/odeint.hpp>

const double LMIN= 0.0;
const double LMAX= 10.0;
const double HPRIME= 0.05;
const double TINF= 200;
const double TMAX= 400;
const double TMIN= 300;
const double DL=0.1;

typedef std::array < double, 2> state_type;

void deriv(const state_type &x, state_type &dxdt, double t);
void write_system(const state_type &x, const double t);

int main(int argc, char**argv)
{
  state_type x={TMIN, 1.3}; //initial conditions

  //-----------------------------funcion,estado,star,fin,paso,
  boost::numeric::odeint::integrate(deriv, x, LMIN, LMAX, DL, write_system);  
}
void deriv(const state_type &x, state_type &dxdt, double t)
{
//x[0]=T y x[1]=z=dT/dx
  dxdt[0]=x[1];
  dxdt[1]=-HPRIME*(TINF-x[0]);
}
void write_system(const state_type &x, const double t)
{
  std::cout << t << "\t" << x[0] << "\t" << x[1] << std::endl;
}
