#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>


using namespace std;
const double ErrMax=1e-7;

double f(double alpha, double r, double l, double t){
  return cos(alpha*t-l*r*sin(t));
}

double IntegralPorSimpson(double alpha, double r, double l, double a, double b, int n){
    double t,h,sum; int i;
    n=2*n; 
    h=(b-a)/n;
    sum=0;
    for(i=0;i<=n;i++){
        t=a+i*h;
        if(i==0 || i==n)
            sum=sum+f(alpha,r,l,t);
        
        else if(i%2==0)
            sum=sum+2*f(alpha,r,l,t);
        
        else
            sum=sum+4*f(alpha,r,l,t);
        
    }
    return sum*h/3;
}

double Bessel(double alpha,double r, double l){
    double a=0,b=M_PI; int n=50;
    return 1.0/M_PI*IntegralPorSimpson(alpha,r,l,a,b,n);
}

double ZerosBi(double alpha,double r, double l, double a, double b){
    double x;
    double m, fa, fm;

    fa=Bessel(alpha,r,a);
    while(b-a >= ErrMax){
        m = (b+a)/2; fm=Bessel(alpha,r,m);
        if(fa*fm>0)
        {a=m; fa=fm;}
        else
        b=m;
        }
        return (a+b)/2;
}

int main(){
    double alpha=0;
    double r=1; 
    double l;
    double zero1=ZerosBi(alpha,r,l,0.01,4);
    double zero2=ZerosBi(alpha,r,l,4,7);
    double zero3=ZerosBi(alpha,r,l,7,10);
    double zero4=ZerosBi(alpha,r,l,10,13);
    double zero5=ZerosBi(alpha,r,l,13,15);
    cout<<"El primer cero está en "<<zero1<<endl;
    cout<<"El segundo cero está en "<<zero2<<endl;
    cout<<"El tercer cero está en "<<zero3<<endl;
    cout<<"El cuarto cero está en "<<zero4<<endl;
    cout<<"El quinto cero está en "<<zero5<<endl;

    ofstream outfile;
    outfile.open("T1_P2_d_zero1.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<Bessel(alpha,r,zero1)<<endl;
    }
    outfile.close();
    outfile.open("T1_P2_d_zero2.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<Bessel(alpha,r,zero2)<<endl;
    }
    outfile.close();
    outfile.open("T1_P2_d_zero3.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<Bessel(alpha,r,zero3)<<endl;
    }
    outfile.close();
    outfile.open("T1_P2_d_zero4.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<Bessel(alpha,r,zero4)<<endl;
    }
    outfile.close();
    outfile.open("T1_P2_d_zero5.dat");
    for(r=0.01;r<=1;r+=0.01){
        outfile <<r<<"\t"<<Bessel(alpha,r,zero5)<<endl;
    }
    outfile.close();
}
