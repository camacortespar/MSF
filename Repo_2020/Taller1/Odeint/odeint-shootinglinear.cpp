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
void shootinglinear(double x0, double xf);

int main(int argc, char**argv)
{
  shootinglinear(TMIN, TMAX); 
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
void shootinglinear(double x0, double xf)
{
  //probar una condicion inicial
  double dev1=-5.0;
  state_type sA={x0, dev1};
 //-----------------------------funcion,estado,star,fin,paso,
 boost::numeric::odeint::integrate(deriv, sA, LMIN, LMAX, DL);
 std::cout << "# final state 1="<<sA[0]<<std::endl;
  //luego otra
  double dev2=-20.0;
  state_type sB={x0, dev2};
 boost::numeric::odeint::integrate(deriv, sB, LMIN, LMAX, DL);
  std::cout << "# final state 2="<<sB[0]<<std::endl;
  //hacer la interpolacion
 double dev = dev1+(dev2-dev1)*(xf-sA[0])/(sB[0]-sA[0]);

  //resolver el sistema con las condiciones iniciales correctas
  state_type s={x0, dev};
  boost::numeric::odeint::integrate(deriv, s, LMIN, LMAX, DL, write_system);
   std::cout << "# dev estimation ="<<dev<<std::endl;
  }
