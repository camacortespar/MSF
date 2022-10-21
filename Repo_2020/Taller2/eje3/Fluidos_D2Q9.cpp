#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

const int Lx=512;
const int Ly=64;

const int Q=9;

const double tau=1.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

const double RHOinicial= 1.0;
const double nu = (tau-0.5)/3;
const double eta = RHOinicial*nu;
const int N = 48;

//----------------------------Clase Lattice Boltzmann --------------------------
class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; //V[0][i]=vi_x, //V[1][i]=vi_y
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q];  //f[ix][iy][i]
  double sigmaxx[Lx][Ly],sigmayy[Lx][Ly],sigmaxy[Lx][Ly];
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool Shownew);
  double Jx(int ix,int iy,bool Shownew);
  double Jy(int ix,int iy,bool Shownew);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Colisione(void);
  void ImponerCampos(double Uventilador);
  void Adveccione(void);
  void Inicie(double rho0,double Ux0,double Uy0);
  void Imprimase(const char *NombreArchivo,double Uventilador);
  void Sigmaxx(double p,double dt);
  void Sigmayy(double p,double dt);
  void Sigmaxy(double p,double dt);
  void ImprimirEsfuerzos(void);
  void Interpolar(double x,double y,double dAx,double dAy,double &Fx,double &Fy);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //cargar los pesos
  w[0]=4.0/9.0; w[1]=w[2]=w[3]=w[4]=1.0/9.0;  w[5]=w[6]=w[7]=w[8]=1.0/36;
  //cargar los vectores velocidad
  V[0][0]=0;  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1;  V[0][4]=0;
  V[1][0]=0;  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;
              V[0][5]=1;  V[0][6]=-1;  V[0][7]=-1;  V[0][8]=1;
              V[1][5]=1;  V[1][6]=1;  V[1][7]=-1;  V[1][8]=-1;
}

double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double suma; int i;
  for (suma=0,i = 0; i < Q; i++) {
    if (UseNew) {
      suma+=fnew[ix][iy][i];
    } else {
      suma+=f[ix][iy][i];
    }
  }
  return suma;
}

double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double suma; int i;
  for (suma=0,i = 0; i < Q; i++) {
    if (UseNew) {
      suma+=V[0][i]*fnew[ix][iy][i];
    } else {
      suma+=V[0][i]*f[ix][iy][i];
    }
  }
  return suma;
}

double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double suma; int i;
  for (suma=0,i = 0; i < Q; i++) {
    if (UseNew) {
      suma+=V[1][i]*fnew[ix][iy][i];
    } else {
      suma+=V[1][i]*f[ix][iy][i];
    }
  }
  return suma;
}

double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i]; double U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Ux0,Uy0;
  for (ix = 0; ix < Lx; ix++) {  // Para cada celda
    for (iy = 0; iy < Ly; iy++) {
      //calcular los campos macroscopicos en la celda
      rho0=rho(ix,iy,false);  Ux0=Jx(ix,iy,false)/rho0;  Uy0=Jy(ix,iy,false)/rho0;
      for (i = 0; i < Q; i++) { //en cada direccion
        fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }
  }
}

void LatticeBoltzmann::ImponerCampos(double Uventilador){
  int ix,iy,i;  double rho0; int ixc=128, iyc=32, R=8; double R2=R*R;
  // me voy por todas las celdas, para imponer si son ventilador u obstaculo
  for (ix = 0; ix < Lx; ix++) {
    for (iy = 0; iy < Ly; iy++) {
      rho0=rho(ix,iy,false);
      if(ix==0){ //ventilador
        for (i = 0; i < Q; i++) {
          fnew[ix][iy][i]=feq(rho0,Uventilador,0,i);
        }
      } else {  //obstaculo
        if ((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2) {
          for (i = 0; i < Q; i++) {
            fnew[ix][iy][i]=feq(rho0,0,0,i);
          }
        }
      }
    }
  }
}

void LatticeBoltzmann::Adveccione(void){
  for (int ix = 0; ix < Lx; ix++) {  // Para cada celda
    for (int iy = 0; iy < Ly; iy++) {
      for (int i = 0; i < Q; i++) { //en cada direccion
        f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i]; // Fronteras periodicas
      }
    }
  }
}

void LatticeBoltzmann::Inicie(double rho0,double Ux0,double Uy0){
  for (int ix = 0; ix < Lx; ix++) {  // Para cada celda
    for (int iy = 0; iy < Ly; iy++) {
      for (int i = 0; i < Q; i++) { //en cada direccion
        f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
      }
    }
  }
}

void LatticeBoltzmann::Imprimase(const char *NombreArchivo,double Uventilador){
  ofstream MiArchivo(NombreArchivo); double rho0,Ux0,Uy0;
  for (int ix = 0; ix < Lx; ix+=4) {  // Para cada celda
    for (int iy = 0; iy < Ly; iy+=4) {
      rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0;  Uy0=Jy(ix,iy,true)/rho0;
      MiArchivo<<ix<<" "<<iy<<" "<<Ux0/Uventilador*4.0<<" "<<Uy0/Uventilador*4.0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Sigmaxx(double p,double dt){
  int ix,iy,i;  double rho0,Ux0,suma;
  for (ix = 0; ix < Lx; ix++) {  // Para cada celda
    for (iy = 0; iy < Ly; iy++) {
      rho0=rho(ix,iy,true);
      for (suma=0,i = 0; i < Q; i++) { //en cada direccion
        Ux0=Jx(ix+V[0][i]*dt,iy+V[1][i]*dt,true)/rho0;
        suma+=w[i]*V[0][i]*Ux0;
      }
      suma*=(3.0/dt);
      sigmaxx[ix][iy]=-p+eta*2.0*suma;
    }
  }
}

void LatticeBoltzmann::Sigmayy(double p,double dt){
  int ix,iy,i;  double rho0,Uy0,suma;
  for (ix = 0; ix < Lx; ix++) {  // Para cada celda
    for (iy = 0; iy < Ly; iy++) {
      rho0=rho(ix,iy,true);
      for (suma=0,i = 0; i < Q; i++) { //en cada direccion
        Uy0=Jy(ix+V[0][i]*dt,iy+V[1][i]*dt,true)/rho0;
        suma+=w[i]*V[1][i]*Uy0;
      }
      suma*=(3.0/dt);
      sigmayy[ix][iy]=-p+eta*2.0*suma;
    }
  }
}

void LatticeBoltzmann::Sigmaxy(double p,double dt){
  int ix,iy,i;  double rho0,Ux0,Uy0,suma;
  for (ix = 0; ix < Lx; ix++) {  // Para cada celda
    for (iy = 0; iy < Ly; iy++) {
      rho0=rho(ix,iy,true);
      for (suma=0,i = 0; i < Q; i++) { //en cada direccion
        Ux0=Jx(ix+V[0][i]*dt,iy+V[1][i]*dt,true)/rho0;  Uy0=Jy(ix+V[0][i]*dt,iy+V[1][i]*dt,true)/rho0;
        suma+=w[i]*(V[1][i]*Ux0+V[0][i]*Uy0);
      }
      suma*=(3.0/dt);
      sigmaxy[ix][iy]=eta*suma;
    }
  }
}

void LatticeBoltzmann::ImprimirEsfuerzos(void){
  int ix,iy,i;
  for (ix = 0; ix < Lx; ix++) {  // Para cada celda
    for (iy = 0; iy < Ly; iy++) {
      cout<<sigmaxx[ix][iy]<<"\t";
    }
    cout<<endl;
  }
  cout<<endl;
}

void LatticeBoltzmann::Interpolar(double x,double y,double dAx,double dAy,double &Fx,double &Fy){
  double u,v; int n=(int)x,m=(int)y; double Deltax, Deltay;
  double sigmaxxp,sigmayyp,sigmaxyp;
  Deltax=1.0; Deltay=1.0;
  // determinar a que celda pertenece el punto que paso
  if (abs(x-n)<=Deltax/2 && abs(y-m)<=Deltay/2) {
    u=(x-n)/Deltax;  v=(y-m)/Deltay;
    // interpolar los cuatro campos en las cuatro celdas adyacentes superiores derechas
    sigmaxxp=sigmaxx[n][m]*(1-u)*(1-v)+sigmaxx[n+1][m]*u*(1-v)+sigmaxx[n][m+1]*(1-u)*v+sigmaxx[n+1][m+1]*u*v;
    sigmayyp=sigmayy[n][m]*(1-u)*(1-v)+sigmayy[n+1][m]*u*(1-v)+sigmayy[n][m+1]*(1-u)*v+sigmayy[n+1][m+1]*u*v;
    sigmaxyp=sigmaxy[n][m]*(1-u)*(1-v)+sigmaxy[n+1][m]*u*(1-v)+sigmaxy[n][m+1]*(1-u)*v+sigmaxy[n+1][m+1]*u*v;
  } else {
    if (abs(x-n)<Deltax/2 && abs(y-m)>Deltay/2) {
      u=(x-n)/Deltax;  v=(y-(m+1))/Deltay;
      // interpolar los cuatro campos en las cuatro celdas adyacentes superiores derechas
      sigmaxxp=sigmaxx[n][m+1]*(1-u)*(1-v)+sigmaxx[n+1][m+1]*u*(1-v)+sigmaxx[n][m+2]*(1-u)*v+sigmaxx[n+1][m+2]*u*v;
      sigmayyp=sigmayy[n][m+1]*(1-u)*(1-v)+sigmayy[n+1][m+1]*u*(1-v)+sigmayy[n][m+2]*(1-u)*v+sigmayy[n+1][m+2]*u*v;
      sigmaxyp=sigmaxy[n][m+1]*(1-u)*(1-v)+sigmaxy[n+1][m+1]*u*(1-v)+sigmaxy[n][m+2]*(1-u)*v+sigmaxy[n+1][m+2]*u*v;
    } else {
      if (abs(x-n)>Deltax/2 && abs(y-m)<Deltay/2) {
        u=(x-(n+1))/Deltax;  v=(y-m)/Deltay;
        // interpolar los cuatro campos en las cuatro celdas adyacentes superiores derechas
        sigmaxxp=sigmaxx[n+1][m]*(1-u)*(1-v)+sigmaxx[n+2][m]*u*(1-v)+sigmaxx[n+1][m+1]*(1-u)*v+sigmaxx[n+2][m+1]*u*v;
        sigmayyp=sigmayy[n+1][m]*(1-u)*(1-v)+sigmayy[n+2][m]*u*(1-v)+sigmayy[n+1][m+1]*(1-u)*v+sigmayy[n+2][m+1]*u*v;
        sigmaxyp=sigmaxy[n+1][m]*(1-u)*(1-v)+sigmaxy[n+2][m]*u*(1-v)+sigmaxy[n+1][m+1]*(1-u)*v+sigmaxy[n+2][m+1]*u*v;
      } else {
        u=(x-(n+1))/Deltax;  v=(y-(m+1))/Deltay;
        // interpolar los cuatro campos en las cuatro celdas adyacentes superiores derechas
        sigmaxxp=sigmaxx[n+1][m+1]*(1-u)*(1-v)+sigmaxx[n+2][m+1]*u*(1-v)+sigmaxx[n+1][m+2]*(1-u)*v+sigmaxx[n+2][m+2]*u*v;
        sigmayyp=sigmayy[n+1][m+1]*(1-u)*(1-v)+sigmayy[n+2][m+1]*u*(1-v)+sigmayy[n+1][m+2]*(1-u)*v+sigmayy[n+2][m+2]*u*v;
        sigmaxyp=sigmaxy[n+1][m+1]*(1-u)*(1-v)+sigmaxy[n+2][m+1]*u*(1-v)+sigmaxy[n+1][m+2]*(1-u)*v+sigmaxy[n+2][m+2]*u*v;
      }
    }
  }
  // calcular la fuerza por unidad de longitud
  Fx=sigmaxxp*dAx+sigmaxyp*dAy;
  Fy=sigmaxyp*dAx+sigmayyp*dAy;
}


//------------------------------------Funciones Globales------------------

int main(void) {
  LatticeBoltzmann Aire;
  int i;
  int t,tmax=100;
  double Uventilador=0.1,dt=1.0,dphi=2*M_PI/N;
  int ixc=128, iyc=32, R=8;
  double xc,yc; //xcircunferencia, ycircunferencia
  double dAx,dAy; // componentes de los diferenciales de superficie
  double Fx,Fy,sumaFx,sumaFy;
  double Re,C_a;


  for (Uventilador = 0; Uventilador < 0.5; Uventilador+=0.01) {
    //Inicie
    Aire.Inicie(RHOinicial,Uventilador,0);
    //Corra
    for (t = 0; t < tmax; t+=dt) {
      Aire.Colisione();
      Aire.ImponerCampos(Uventilador);
      Aire.Adveccione();
      Aire.Sigmaxx(RHOinicial/3,dt);
      Aire.Sigmayy(RHOinicial/3,dt);
      Aire.Sigmaxy(RHOinicial/3,dt);
      for (sumaFx=0,sumaFy=0,i = 1; i <= N; i++) {
        xc=ixc+R*cos(dphi*i);  yc=iyc+R*sin(dphi*i);
        dAx=R*dphi*cos(dphi*i);   dAy=R*dphi*sin(dphi*i);
        Aire.Interpolar(xc,yc,dAx,dAy,Fx,Fy);
        sumaFx+=Fx;  sumaFy+=Fy;
      }
      //cout<<sumaFx<<"\t"<<sumaFy<<endl;
    }
    Re=2*R*Uventilador/nu; C_a=(2*sumaFx)/(RHOinicial*2*R*Uventilador);
    cout<<Re<<"\t"<<C_a<<endl;
  }
  Aire.Imprimase("Aire.dat",Uventilador);


  return 0;
}
