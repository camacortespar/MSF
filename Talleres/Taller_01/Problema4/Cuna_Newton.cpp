#include <iostream>
#include <cmath>

using namespace std;

void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);

// Constantes globales en CGS
const int N = 3; // # péndulos
const double g = 980;
const double K = 10.0e10;

// Constantes PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

// Clases
class Cuerpo;
class Colisionador;

class Cuerpo{
private:
  double theta, omega, tau;
  double m, R, l, x0, I;
public:
  void Inicie(double theta0, double omega0, double m0, double R0, double l0, double x00);
  void BorreTorque(void){tau = 0;};
  void SumeTorque(double tau0){tau += tau0;};
  void Mueva_theta(double dt, double coeficiente);
  void Mueva_omega(double dt, double coeficiente);
  void Dibujese(void);
  double Getx(void){return x0 + l*sin(theta);};
  double Gety(void){return -l*cos(theta);};
  double Gettau(void){return tau;};
  friend class Colisionador;
};

void Cuerpo::Inicie(double theta0, double omega0, double m0, double R0, double l0, double x00)
{
  theta = theta0;
  omega = omega0;
  l = l0;
  x0 = x00;
  m=m0;
  R=R0;
  I=m*l*l;
}

void Cuerpo::Mueva_theta(double dt, double coeficiente)
{
  theta += omega*(dt*coeficiente);
}

void Cuerpo::Mueva_omega(double dt, double coeficiente)
{
  omega += tau*(dt*coeficiente/I);
}

void Cuerpo::Dibujese(void)
{
  // Dibuja bola
  cout << " , " << Getx() << "+" << R << "*cos(t)," << Gety() << "+" << R << "*sin(t)";
  // Dibuja hilo
  cout << " , " << x0 << "+" << l/7 << "*t*sin("<<theta<<"),-" << l/7 << "*t*cos("<<theta<<")";
}

class Colisionador{
private:
public:
  void CalculeTorques(Cuerpo * Pendulo);
  void CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2nn);
};

void Colisionador::CalculeTorques(Cuerpo * Pendulo)
{
  int i;
  double tau0;
  // Actualiza Torques
  for (i = 0; i < N; ++i) {
    Pendulo[i].BorreTorque();
    tau0 = -Pendulo[i].l*Pendulo[i].m*g*sin(Pendulo[i].theta);
    Pendulo[i].SumeTorque(tau0);
  }
  // Interacción entre pendulos
  for (i = N - 1; i > 0; i--) {
    CalculeTorqueEntre(Pendulo[i], Pendulo[i - 1]);
  }
}

void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2)
{
  double s = (Pendulo2.Getx() + Pendulo2.R) - (Pendulo1.Getx()- Pendulo1.R);
  double F = 0;
  if (s > 0) {
    F = K*pow(s, 1.5); // Fuerza de Hertz
  }
  Pendulo1.SumeTorque(F*Pendulo1.l);
  Pendulo2.SumeTorque(-F*Pendulo2.l);
}

// ----------------------------------------------------
int main()
{
  Cuerpo Pendulo[N];
  Colisionador Hertz;
  double m0 = 100, l0 = 12, R0 = 1.5 ;
  double T = 2*M_PI*sqrt(l0/g);
  double t, tmax = 1.0*T, dt = 0.000001;
  double tdibujo, tcuadro = T/100;
  int i;

  // Inicializar pendulos
  //---------------(theta0,omega0,m0,R0,l0,x0)
  Pendulo[0].Inicie(-0.2617993878, 0, m0, R0, l0, 0);
  for (i = 1; i < N; ++i) {
    Pendulo[i].Inicie(0, 0, m0, R0, l0, 2*R0*i);
  }

  // Dibujos
  //InicieAnimacion();

   for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
    //Dibujar
     // if(tdibujo>tcuadro){
     //   InicieCuadro();
     //   for(i=0;i<N;i++) Pendulo[i].Dibujese();
     //   TermineCuadro();
     //   tdibujo=0;
     // }
     
    
    //cout<<Pendulo[1].Getx()<<" "<<Pendulo[1].Gety()<<endl;
    // Mover por PEFRL
     for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);
     Hertz.CalculeTorques(Pendulo);
     for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
     for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
     Hertz.CalculeTorques(Pendulo);
     for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
     for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Coeficiente2);
     Hertz.CalculeTorques(Pendulo);
     for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
     for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
     Hertz.CalculeTorques(Pendulo);
     for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
     for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);

     //P4a, Torque del pendulo de la mitad
     cout << t << " " << Pendulo[1].Gettau() << endl;
   }

   return 0;
}
// ----------------------------------------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'TresPendulos.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:15]"<<endl;
  cout<<"set yrange[-15:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}


