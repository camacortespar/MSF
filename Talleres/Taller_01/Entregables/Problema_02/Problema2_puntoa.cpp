#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;

//Problem 2 - a .
//Solving the dif eq, using 4th order Runge Kutta.

//We write the equation as two coupled equations

double f1(double r, double l, double R1, double R2){
  return -R1/r-(l*l*R2); 
}

double f2(double r, double l, double R1, double R2){
    return R1;
}

void RungeKuttaStep4(double & r0, double & l0, double & R10, double & R20, double dr){  
    double dR11, dR21, dR31, dR41; 
    double dR12, dR22, dR32, dR42;
    dR11=f1(r0,l0,R10,R20)*dr;                       dR12=f2(r0,l0,R10,R20)*dr;
    dR21=f1(r0+dr/2,l0,R10+dR11/2,R20+dR12/2)*dr;    dR22=f2(r0+dr/2,l0,R10+dR11/2,R20+dR12/2)*dr;
    dR31=f1(r0+dr/2,l0,R10+dR21/2,R20+dR22/2)*dr;    dR32=f2(r0+dr/2,l0,R10+dR21/2,R20+dR22/2)*dr;
    dR41=f1(r0+dr,l0,R10+dR31,R20+dR32)*dr;          dR42=f2(r0+dr,l0,R10+dR31,R20+dR32)*dr;

    R10+=(dR11+2*(dR21+dR31)+dR41)/6;             R20+=(dR12+2*(dR22+dR32)+dR42)/6;
    r0+=dr;
}

int main(){
    double r,l,R1,R2; 
    double dr=0.01;
    ofstream outfile;
    outfile.open("T1_P2_av2.dat");
    for(r=0.01,l=1,R1=0,R2=1;r<10;r+=0.01){
        outfile <<r<<" "<<R2<<endl;
        RungeKuttaStep4(r,l,R1,R2,dr);
    }
    outfile.close(); 
return 0; 
}