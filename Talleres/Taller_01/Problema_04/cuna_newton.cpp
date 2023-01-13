//Ejercicio 4.a

#include <iostream>
#include <cmath>

using namespace std;

//Constantes globales (todo en CGS)
const int N = 3;        //numero de pendulos
const double g = 980; 
const double K = 10e12;  //cte fuerza de Hertz

//Constantes PEFRL
const double zeta = 0.1786178958448091e00;
const double lambda = -0.2123418310626054e0;
const double chi = -0.6626458266981849e-1;
const double coef1 = (1-2*lambda)/2;
const double coef2 = 1-2*(chi+zeta);

//Clases
class cuerpo;
class colisionador;

//----------CUERPO----------
//Descripcion de los pendulos
class cuerpo{  
private:  
  double theta, omega, tau, I;  //cantidades angulares
  double m, R, l, x0;           //x0 es distancia horizontal medida desde el origen
public:
  void inicie(double theta0, double omega0, double m0, double R0, double l0, double xi);
  void borretorque(void){tau = 0;};
  void addtorque(double tau0){tau += tau0;};
  void mueva_theta(double dt, double coef);
  void mueva_omega(double dt, double coef);
  void dibujar(void);                          //animacion
  double getx(void){return x0+l*sin(theta);};  //informacion de los pendulos para dibujar 
  double gety(void){return -l*cos(theta);};
  double gettau(void){return tau;};
  //para que todo lo del colisionador pueda meterle la mano a los datos privados de cuerpo
  friend class colisionador;
};

//funciones
void cuerpo::inicie(double theta0, double omega0, double m0, double R0, double l0, double xi){
  theta = theta0; omega = omega0; m = m0; R = R0; l = l0; x0 = xi; I = m*pow(l,2);
}

void cuerpo::mueva_theta(double dt, double coef){
  theta += omega*(coef*dt);
}

void cuerpo::mueva_omega(double dt, double coef){
  omega += tau*(coef*dt/I);
}

//animacion
void cuerpo::dibujar(void){
  cout<<" , "<<getx()<<"+"<<R<<"*cos(t),"<<gety()<<"+"<<R<<"*sin(t)";              //esferas
  cout<<" , "<<x0<<"+"<<l/7<<"*t*sin("<<theta<<"),-"<<l/7<<"*t*cos("<<theta<<")";  //cuerdas
}

//----------COLISIONADOR----------
class colisionador{  
private:
public:
  void calculetorques(cuerpo * pend);
  void entretorque(cuerpo & pend1, cuerpo & pend2);
};

//funciones
void colisionador::calculetorques(cuerpo * pend){
  int i, j; double tau0;
  for(i = 0; i < N; i++){
    pend[i].borretorque();                             //borra torques 
    tau0 = -pend[i].l*pend[i].m*g*sin(pend[i].theta);  //torque que hace la fuerza de gravedad
    pend[i].addtorque(tau0);                           //adicione torques
  }
  //calcule torque entre todas las parejas
  for(i = N-1; i > 0; i--){
    entretorque(pend[i], pend[i-1]);  //pendulo de la derecha interactua con el de la izquierda
  }
}

void colisionador::entretorque(cuerpo & pend1, cuerpo & pend2){
  double s = (pend2.getx()+pend2.R)-(pend1.getx()-pend1.R), f = 0;
  if(s > 0){
    f = K*pow(s,1.5);  //fuerza de Hertz
  }
  pend1.addtorque(f*pend1.l); pend2.addtorque(-f*pend2.l);  //torque que hace la fuerza de Hertz
}

//----------ANIMACION----------
void animacion(void){
  //cout<<"set terminal gif animate"<<endl;  //realizar un gif
  //cout<<"set output 'pends.gif'"<<end;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-5:11]"<<endl;           //definir rangos
  cout<<"set yrange[-15:0]"<<endl;
  cout<<"set size ratio -1"<<endl;           //proporcion 1:1 entre las unidades de x e y
  cout<<"set parametric"<<endl;              //se realiza una grafica parametrica
  cout<<"set trange [0:7]"<<endl;            //rango de la animacion
  cout<<"set isosamples 12"<<endl;           //hace 12 divisiones dentro del rango  
}

void iniciecuadro(void){
  cout<<"plot 0,0 ";
}

void terminecuadro(void){
  cout<<endl;
}

//---------------//
//Tres pendulos interactuando
int main(void){
  cuerpo pend[N];
  colisionador hertz;
  double m0 = 100, l0 = 12, R0 = 1.5, the0 = -0.2617;  //the10 es el angulo inicial del pend[0]
  double T=2*M_PI*sqrt(l0/g);                          //periodo de oscilacion de un pendulo
  double t, tmax = (0.5)*T, dt = 0.000001;  
  double tdibujo, tcuadro = T/5000; 
  int i;
  
  //inicializacion
  //------------(theta0, omega0, m0, R0, l0, xi)
  pend[0].inicie(the0,        0, m0, R0, l0, 0);  //primer pendulo separado
  for(i = 1; i < N; i++){
    pend[i].inicie(0, 0, m0, R0, l0, 2*R0*i);
  }

  //animacion();
  for(t = 0, tdibujo = 0; t < tmax; t += dt, tdibujo += dt){
    /*if(tdibujo > tcuadro){
      iniciecuadro();
      for(i = 0; i < N; i++){
	pend[i].dibujar();
      }
      terminecuadro();
      tdibujo = 0;
      }*/
    
    cout<<t<<" "<<pend[1].gettau()<<endl;  //datos del pend intermedio
    
    //movimiento por PEFRL
    for(i = 0; i < N; i++){pend[i].mueva_theta(dt, zeta);}
    hertz.calculetorques(pend);
    for(i = 0; i < N; i++){pend[i].mueva_omega(dt, coef1);}
    for(i = 0; i < N; i++){pend[i].mueva_theta(dt, chi);}
    hertz.calculetorques(pend);
    for(i = 0; i < N; i++){pend[i].mueva_omega(dt, lambda);}
    for(i = 0; i < N; i++){pend[i].mueva_theta(dt, coef2);}
    hertz.calculetorques(pend);
    for(i = 0; i < N; i++){pend[i].mueva_omega(dt, lambda);}
    for(i = 0; i < N; i++){pend[i].mueva_theta(dt, chi);}
    hertz.calculetorques(pend);
    for(i = 0; i < N; i++){pend[i].mueva_omega(dt, coef1);}
    for(i = 0; i < N; i++){pend[i].mueva_theta(dt, zeta);}
  }   

  return 0;
}

  
