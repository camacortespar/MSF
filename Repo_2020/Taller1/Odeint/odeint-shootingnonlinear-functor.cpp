//g++ -std=c++11 odeint-shootingnonlinear.cpp
#include <iostream>
#include <array>
#include <boost/numeric/odeint.hpp>

const double LMIN= 0.0;
const double LMAX= 10.0;
const double HPRIME= 0.05;
const double SIGMA=2.7e-10;
const double TINF= 200;
const double TMAX= 400;
const double TMIN= 300;
const double DL=0.1;

typedef std::array < double, 2> state_type;

void deriv(const state_type &x, state_type &dxdt, double t);
void write_system(const state_type &x, const double t);
void shootinglinear(double x0, double xf);
void shooting_non_linear(double x0, double xf);
template <typename F>
double bisection(double xl0, double xu0, double eps, F f);


int main(int argc, char**argv)
{
  shootinglinear(TMIN, TMAX);
  std::cout <<"\n\n";
  shooting_non_linear(TMIN, TMAX);
}
void deriv(const state_type &x, state_type &dxdt, double t)
{
//x[0]=T y x[1]=z=dT/dx
  dxdt[0]=x[1];
  dxdt[1]=-HPRIME*(TINF-x[0])-SIGMA*(std::pow(TINF,4)-std::pow(x[0],4));
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

struct Functor{
  double x0=0;
  double xf=0;
  double operator() (double z){
    state_type s={x0, z};
    boost::numeric::odeint::integrate(deriv, s, LMIN, LMAX, DL);
    return s[0]-xf;
  };
};



  void shooting_non_linear(double x0, double xf)
  {
  //biseccion con la funcion que compare xf[0]-xf
  //usando una funcion lambda [](){}
  //double dev=bisection(-5.0, -20.0, 1.0e-4,[x0,xf](double z){
    //					     state_type s={x0, z};
    //					     boost::numeric::odeint::integrate(deriv, s, LMIN, LMAX, DL);
    //					     return s[0]-xf;
    //					   });
    Functor fun;
    fun.x0=x0;
    fun.xf=xf;
    double dev=bisection(-5.0, -20.0, 1.0e-4, fun);
    
    state_type s={x0, dev};
  boost::numeric::odeint::integrate(deriv, s, LMIN, LMAX, DL, write_system);
}
template <typename F>
double bisection(double xl0, double xu0, double eps, F f)
{
  int iter = 1;
  double mult=f(xl0)*f(xu0);
  double xl=xl0;
  double xu=xu0;
  
  while (mult < 0){
    double xr = (xl + xu)/2;
    //std::cout << iter << "\t"<<xl<< <<xu<< <<xr<< fun(xr)<<"\n";
    if(std::fabs(f(xr))<=eps){
      return xr;
    }
    if(f(xr)*f(xu)<0){
      xl=xr;
    } else {
      xu=xr;
    }
    mult=f(xu)*f(xl);
    iter++;
  }
  return 0.5*(xl+xu);
}
