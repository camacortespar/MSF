//planetary_data.h
#include"master.h"

const double G_meters=6.67408e-11;
const double G=G_meters*std::pow(10,-9);

class sun {
 private:

 public:
  std::string name="SUN";
  double mass=1.989e30;
  double mu=1.32712e11;
  double radius=695700.0;
  double g0=274.0;
  double J2=-6.84e-7;
  std::vector<double> atmrotvector{0,0,72.9211e-6}; //dato de la tierra
  std::vector<double> zs{63.096,251.189,1000}; //km dato de la tierra
  std::vector<double> rhos{2.059e5,5.909e-2,3.561e-6};// # kg/km³ dato de la tierra
  double deorbit_altitude=1000;
};

class earth {
 private:

 public:
  std::string name="EARTH";
  double mass=5.972e24;
  double mu=5.972e24*G;
  double radius=6378.0;
  double J2=1.082635854e-3;
  double g0=9.81;
  std::vector<double> atmrotvector{0,0,72.9211e-6};
  std::vector<double> zs{63.096,251.189,1000}; //km
  std::vector<double> rhos{2.059e5,5.909e-2,3.561e-6};// # kg/km³
  double deorbit_altitude=100;
  //'spice_file':A+'de432s.bsp'
};

class perturbations {
 private:

 public:
  bool J2=false;
  bool aero=false;
  bool srp=false;
  //n_bodies=[];
  double Cd=0;
  double A=0;
  double thrust=0;
  double thrust_direction=0;
  double isp=0;
  double CR=0;
  double A_srp=0;
  double C20=0;
  double mu=0;
};
class StopC  //definir las condiciones para que apague el motor
{
private:

public:
  double max_alt=2000;
  double min_alt=0;
  bool no=false;
};

