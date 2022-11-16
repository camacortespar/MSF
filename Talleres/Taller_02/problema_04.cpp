/*
  Taller 02 - Ejercicio 04:
  Fuerza de Magnus
*/

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//-----Constantes Globales-----//
const int Lx = 512;                 //tamaño de la simulacion
const int Ly = 64;

const int Q = 9;                    //numero de direcciones

const double tau = 1.5;             //valores en la funcion de colision
const double Utau = 1.0/tau;
const double UmUtau = 1-Utau;
//Propiedades del fluido
const double rho_fluido = 1.0;
const double nu_fluido = (tau-0.5)/3.0;           //viscosidad cinematica
const double eta_fluido = rho_fluido*nu_fluido;   //viscosidad dinamica 

//-----Clase LatticeBoltzmann------
class LatticeBoltzmann{
private:
  double w[Q];                                                //pesos por direccion
  int Vx[Q], Vy[Q];                                           //vectores de velocidad
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q];                       //funciones de distribucion
  double sigmaxx[Lx][Ly], sigmayy[Lx][Ly], sigmaxy[Lx][Ly];   //tensor de esfuerzos
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double feq(double rho0, double Ux0, double Uy0, int i);
  void Start(double rho0, double Ux0, double Uy0);
  void Collision(void);
  void ImposeFields(double Ufan, double omega);
  void Advection(void);
  void Sigmaxx(double p0, double eta0, double dt);
  void Sigmayy(double p0, double eta0, double dt);
  void Sigmaxy(double p0, double eta0, double dt);
  void Interpol(double x, double y, double dAx, double dAy, double & Fx, double & Fy);
  void Print(const char * NameFile, double Ufan);
};

//Constructor
LatticeBoltzmann::LatticeBoltzmann(void){   
   //Pesos
  w[0] = 4.0/9.0;  w[1] = w[2] = w[3] = w[4] = 1.0/9.0;  w[5] = w[6] = w[7] = w[8] = 1.0/36.0;
  //Vectores de velocidad
  Vx[0] = 0;  Vx[1] = 1;  Vx[2] = 0;  Vx[3] = -1; Vx[4] = 0;
  Vy[0] = 0;  Vy[1] = 0;  Vy[2] = 1;  Vy[3] = 0;  Vy[4] = -1;
  Vx[5] = 1;  Vx[6] = -1; Vx[7] = -1; Vx[8] = 1;
  Vy[5] = 1;  Vy[6] = 1;  Vy[7] = -1; Vy[8] = -1;
}
//Rho
double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  double sum; int i;
  for(sum=0, i=0; i<Q; i++){
    if(UseNew) sum += fnew[ix][iy][i]; else sum += f[ix][iy][i];
  }
  return sum;
} 
//Componentes del vector J
double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  double sum; int i;
  for(sum=0, i=0; i<Q; i++){
    if(UseNew) sum += Vx[i]*fnew[ix][iy][i]; else sum += Vx[i]*f[ix][iy][i];
  }
  return sum;
}
double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  double sum; int i;
  for(sum=0, i=0; i<Q; i++){
    if(UseNew) sum += Vy[i]*fnew[ix][iy][i]; else sum += Vy[i]*f[ix][iy][i];
  }
  return sum;
}
//Funcion equilibrio
double  LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i], U2 = Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
  } 
//Inicie
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix, iy, i;
  for(ix=0; ix<Lx; ix++)    //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){   //en cada direccion
        f[ix][iy][i] = feq(rho0,Ux0,Uy0,i);
      }
}  
//Colision
void LatticeBoltzmann::Collision(void){
  int ix, iy, i; double rho0, Ux0, Uy0;
  for(ix=0; ix<Lx; ix++)      //para cada celda
    for(iy=0; iy<Ly; iy++){
      //Calcule los campos macroscopicos en la celda
      rho0 = rho(ix,iy,false); Ux0 = Jx(ix,iy,false)/rho0; Uy0 = Jy(ix,iy,false)/rho0;
      for(i=0; i<Q; i++){     //para cada vector de velocidad
        fnew[ix][iy][i] = UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }  
}
//Imponer campos, en este caso se impone geometría de ventilador + cilindro
void LatticeBoltzmann::ImposeFields(double Ufan, double omega){
  int i, ix, iy;
  double rho0; int ixc = 128, iyc = 32, R = 8; double R2 = R*R;   //ctes que definen cilindro
  for(ix=0; ix<Lx ; ix++)      //para cada celda
    for(iy=0; iy<Ly; iy++){
      rho0 = rho(ix, iy, false);
      //Ventilador
      if(ix==0)
        for(i=0; i<Q; i++){fnew[ix][iy][i] = feq(rho0,Ufan,0,i);}
      //Cilindro rotatorio
      else if(pow((ix-ixc),2)+pow((iy-iyc),2) <=R2)
        for(i=0; i<Q; i++){fnew[ix][iy][i] = feq(rho0,-omega*(iy-iyc),omega*(ix-ixc),i);}
      //Hay que agregar un punto de perturbacion para romper la isotropia
      //else if(ix==ixc && iy==iyc+R+1)
      //  for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}
    }
}
//Adveccion
void LatticeBoltzmann::Advection(void){
  int ix, iy, i;
  for(ix=0; ix<Lx; ix++)      //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){     //en cada direccion
        f[(ix+Vx[i]+Lx)%Lx][(iy+Vy[i]+Ly)%Lx][i]= fnew[ix][iy][i];       //fronteras periodicas
      }
}
//Tensor de esfuerzos
void LatticeBoltzmann::Sigmaxx(double p0, double eta0, double dt){
  int ix, iy, i; double rho0, Ux0, sum;
  for(ix=0; ix<Lx; ix++){               //para cada celda
    for(iy=0; iy<Ly; iy++){
      rho0 = rho(ix,iy,true);
      for(sum=0, i=0; i<Q; i++){      //en cada direccion
        Ux0 = Jx(ix+Vx[i]*dt,iy+Vy[i]*dt,true)/rho0;        //velocidad en el siguiente paso
        sum += w[i]*Vx[i]*Ux0;                              //sumatoria en la derivada parcial
      }
      sum *= (3.0/dt);                                      //derivada parcial
      sigmaxx[ix][iy] = -p0+eta0*2.0*sum;                   //sigma xx
    }
  }
}
void LatticeBoltzmann::Sigmayy(double p0, double eta0, double dt){
  int ix, iy, i; double rho0, Uy0, sum;
  for(ix=0; ix<Lx; ix++){               //para cada celda
    for(iy=0; iy<Ly; iy++){
      rho0 = rho(ix,iy,true);
      for(sum=0, i=0; i<Q; i++){      //en cada direccion
        Uy0 = Jy(ix+Vx[i]*dt,iy+Vy[i]*dt,true)/rho0;        //velocidad en el siguiente paso
        sum += w[i]*Vy[i]*Uy0;                              //sumatoria en derivada parcial
      }
      sum *= (3.0/dt);                                      //derivada parcial
      sigmayy[ix][iy] = -p0+eta0*2.0*sum;                   //sigma yy
    }
  }
}
void LatticeBoltzmann::Sigmaxy(double p0, double eta0, double dt){
  int ix, iy, i; double rho0, Ux0, Uy0, sum;
  for(ix=0; ix<Lx; ix++){               //para cada celda
    for(iy=0; iy<Ly; iy++){
      rho0 = rho(ix,iy,true);
      for(sum=0, i=0; i<Q; i++){     //en cada direccion
        Ux0 = Jx(ix+Vx[i]*dt,iy+Vy[i]*dt,true)/rho0; Uy0=Jy(ix+Vx[i]*dt,iy+Vy[i]*dt,true)/rho0;
        sum += w[i]*(Vy[i]*Ux0+Vx[i]*Uy0);
      }
      sum *= (3.0/dt);
      sigmaxy[ix][iy] = eta0*sum;          //sigma xy/yx
    }
  }
}
//Interpolacion
void LatticeBoltzmann::Interpol(double x, double y, double dAx, double dAy, double & fx_aux, double & fy_aux){
  double u, v; int nx = int(x), ny = int(y); double Dx, Dy;
  double sigmaxx_p, sigmayy_p, sigmaxy_p;
  Dx = 1.0; Dy = 1.0;
  //Determine la celda correspondiente al punto P
  if(abs(x-nx)<=Dx/2 && abs(y-ny)<=Dy/2){
    u = (x-nx)/Dx; v = (y-ny)/Dy;
    //Interpolacion bilinieal
    sigmaxx_p = sigmaxx[nx][ny]*(1-u)*(1-v)+sigmaxx[nx+1][ny]*u*(1-v)+sigmaxx[nx][ny+1]*(1-u)*v+sigmaxx[nx+1][ny+1]*u*v;
    sigmayy_p = sigmayy[nx][ny]*(1-u)*(1-v)+sigmayy[nx+1][ny]*u*(1-v)+sigmayy[nx][ny+1]*(1-u)*v+sigmayy[nx+1][ny+1]*u*v;
    sigmaxy_p = sigmaxy[nx][ny]*(1-u)*(1-v)+sigmaxy[nx+1][ny]*u*(1-v)+sigmaxy[nx][ny+1]*(1-u)*v+sigmaxy[nx+1][ny+1]*u*v;
  }
  else{
    if(abs(x-nx)<Dx/2 && abs(y-ny)>Dy/2){
      u = (x-nx)/Dx; v = (y-(ny+1))/Dy;
      //
      sigmaxx_p = sigmaxx[nx][ny+1]*(1-u)*(1-v)+sigmaxx[nx+1][ny+1]*u*(1-v)+sigmaxx[nx][ny+2]*(1-u)*v+sigmaxx[nx+1][ny+2]*u*v;
      sigmayy_p = sigmayy[nx][ny+1]*(1-u)*(1-v)+sigmayy[nx+1][ny+1]*u*(1-v)+sigmayy[nx][ny+2]*(1-u)*v+sigmayy[nx+1][ny+2]*u*v;
      sigmaxy_p = sigmaxy[nx][ny+1]*(1-u)*(1-v)+sigmaxy[nx+1][ny+1]*u*(1-v)+sigmaxy[nx][ny+2]*(1-u)*v+sigmaxy[nx+1][ny+2]*u*v;
    } 
    else{
      if(abs(x-nx)>Dx/2 && abs(y-ny)<Dy/2){
        u = (x-(nx+1))/Dx;  v = (y-ny)/Dy;
        //
        sigmaxx_p = sigmaxx[nx+1][ny]*(1-u)*(1-v)+sigmaxx[nx+2][ny]*u*(1-v)+sigmaxx[nx+1][ny+1]*(1-u)*v+sigmaxx[nx+2][ny+1]*u*v;
        sigmayy_p = sigmayy[nx+1][ny]*(1-u)*(1-v)+sigmayy[nx+2][ny]*u*(1-v)+sigmayy[nx+1][ny+1]*(1-u)*v+sigmayy[nx+2][ny+1]*u*v;
        sigmaxy_p = sigmaxy[nx+1][ny]*(1-u)*(1-v)+sigmaxy[nx+2][ny]*u*(1-v)+sigmaxy[nx+1][ny+1]*(1-u)*v+sigmaxy[nx+2][ny+1]*u*v;
      } 
      else{
        u = (x-(nx+1))/Dx;  v = (y-(ny+1))/Dy;
        //
        sigmaxx_p = sigmaxx[nx+1][ny+1]*(1-u)*(1-v)+sigmaxx[nx+2][ny+1]*u*(1-v)+sigmaxx[nx+1][ny+2]*(1-u)*v+sigmaxx[nx+2][ny+2]*u*v;
        sigmayy_p = sigmayy[nx+1][ny+1]*(1-u)*(1-v)+sigmayy[nx+2][ny+1]*u*(1-v)+sigmayy[nx+1][ny+2]*(1-u)*v+sigmayy[nx+2][ny+2]*u*v;
        sigmaxy_p = sigmaxy[nx+1][ny+1]*(1-u)*(1-v)+sigmaxy[nx+2][ny+1]*u*(1-v)+sigmaxy[nx+1][ny+2]*(1-u)*v+sigmaxy[nx+2][ny+2]*u*v;
      }
    }
  }
  //Calcule fuerza
  fx_aux = sigmaxx_p*dAx+sigmaxy_p*dAy;
  fy_aux = sigmaxy_p*dAx+sigmayy_p*dAy;
}
//Print
void LatticeBoltzmann::Print(const char * NameFile, double Ufan){
  ofstream MyFile(NameFile); double rho0, Ux0, Uy0; int ix, iy;
  for(ix = 0; ix < Lx; ix += 4){
    for(iy = 0; iy < Ly; iy += 4){
      rho0 = rho(ix,iy,true); Ux0 = Jx(ix,iy,true)/rho0; Uy0 = Jy(ix,iy,true)/rho0;
      MyFile<<ix<<" "<<iy<<" "<<Ux0/Ufan*4.0<<" "<<Uy0/Ufan*4.0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//-----Programa Principal-----
int main(void){
  LatticeBoltzmann Air;
  int i;
  int t, tmax = 100; double dt = 1.0;   //tiempos de simulacion
  double Ufan0 = 0.1;                   //propiedades de ventilador + cilindro
  double omega = 2*M_PI/1000;
  int ixc = 128, iyc = 32, R = 8;
  double xc, yc; 
  int N = 24;                           //division en N elementos planos del cilindro
  double dAx, dAy, dphi = 2*M_PI/N;     //componentes del dA
  double fx_aux, fy_aux, Fx, Fy;
  double Re, Ca, A = 2*R, Cm = 1.0; 
  
  for(omega = 0.01; omega <= 0.3 ; omega += 0.01){
    //Inicie
    Air.Start(rho_fluido,Ufan0,0);
    //Evolucione
    for(t = 0; t < tmax; t += dt){
      Air.Collision();
      Air.ImposeFields(Ufan0, omega);
      Air.Advection();
      Air.Sigmaxx(rho_fluido/3.0,eta_fluido,dt);
      Air.Sigmayy(rho_fluido/3.0,eta_fluido,dt);
      Air.Sigmaxy(rho_fluido/3.0,eta_fluido,dt);
      if(t=tmax-dt){
        for(Fx = 0, Fy = 0, i = 1; i <= N; i++){
            xc = ixc+R*cos(dphi*i);  yc = iyc+R*sin(dphi*i);        //encuentre el punto central de dA
            dAx = R*dphi*cos(dphi*i);   dAy = R*dphi*sin(dphi*i);   //componentes de dA
            Air.Interpol(xc, yc, dAx, dAy, fx_aux, fy_aux);
            Fx += fx_aux; Fy += fy_aux;
        }
      }
    }
    //Numero de Reynolds y coeficiente de arrastre
    Re = A*Ufan0/nu_fluido; Ca = (2*Fx)/(rho_fluido*A*Ufan0*Ufan0);
    //cout<<Re<<"\t"<<Ca<<endl;
    cout<<omega<<"\t"<<-Fy<<endl;
  }
  //Air.Print("Air.dat", Ufan0);
  return 0;
}