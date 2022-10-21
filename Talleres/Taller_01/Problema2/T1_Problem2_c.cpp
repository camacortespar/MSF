#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;

//Problem 2 - c .


const double ErrMax=1e-7;

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
   
}

double R(double rf,double lf){
    double r,l,R1,R2,result=0; 
    double dr=0.01;
    for(l=0.1;l<=15;l+=0.01){
        for(r=0.01,R1=0,R2=1;r<=10;r+=0.01){
            
            if(r>rf-0.01 && r<rf+0.01){
                if(l>lf-0.01 && l<lf+0.01)
                    result=R2;
            }
            RungeKuttaStep4(r,l,R1,R2,dr);
        }
    }
    
return result; 
}

double ZerosBi(double a, double b){
    double m, fa, fm;
    //Calculated for a function with r=1
    fa=R(1,a);
    while(b-a >= ErrMax){
        m = (b+a)/2; fm=R(1,m);
        if(fa*fm>0)
        {a=m; fa=fm;}
        else
        b=m;
        }
    return (a+b)/2;
}

int main(){
    double zero1=ZerosBi(0.5,4);
    double zero2=ZerosBi(4,7);
    double zero3=ZerosBi(7,10);
    double zero4=ZerosBi(10,13);
    double zero5=ZerosBi(13,15);
    cout<<"El primer cero está en "<<zero1<<endl;
    cout<<"El segundo cero está en "<<zero2<<endl;
    cout<<"El tercer cero está en "<<zero3<<endl;
    cout<<"El cuarto cero está en "<<zero4<<endl;
    cout<<"El quinto cero está en "<<zero5<<endl;
    ofstream outfile;
    outfile.open("T1_P2_c_zero1.dat");
    double r;
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<R(r,zero1)<<endl;
    }
    outfile.close(); 
    outfile.open("T1_P2_c_zero2.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<R(r,zero2)<<endl;
    }
    outfile.close(); 
    outfile.open("T1_P2_c_zero3.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<R(r,zero3)<<endl;
    }
    outfile.close(); 
    outfile.open("T1_P2_c_zero4.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<R(r,zero4)<<endl;
    }
    outfile.close(); 
    outfile.open("T1_P2_c_zero5.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<R(r,zero5)<<endl;
    }
    outfile.close();

}