//tools.h
#pragma once
#include "master.h"

const double d2r=M_PI/180;
const double r2d=180/M_PI;

//---------------------declarar funciones-----------------------

template <typename T>
void Plot_orbit_gnuplot(std::vector<OrbitPropagator> &OP,const T &cb, std::string title, bool save);
template <typename T>
std::vector<double> coes2rv(std::vector<double> &coes, bool deg, const T& cb);
double ecc_anomaly(double ta, double e, std::string method);
double true_anomaly(double E,double e);
template <typename T>
std::vector<double> tlecoes(std::string file, const T &cb);
template <typename T>
std::vector <double> rv2coes(std::vector <double> state,const T &CB,bool deg);
template <typename T>
void Plot_orbit(std::vector<OrbitPropagator>&OP,T &cb,bool show,bool save,std:: string title,std::vector<std::string> &labels);
template <typename T>
double calc_atm_density(double z, T &CB);
template <typename T>
std::vector <double> find_rho_z(double z, T &CB);
template <typename T>
void Plot_coes(std::vector<OrbitPropagator> &OP,T &CB,bool show, bool save,std::vector<std::string>labels, bool hours, bool days);
template <typename T>
void  calculate_coes(std::vector <OrbitPropagator> &OP,T &CB);

//--------------------------implementar funciones----------------

template <typename T>
void Plot_orbit_gnuplot(std::vector<OrbitPropagator> &OP,const T &cb, std::string title, bool save)
{
 double r=cb.radius;
 int t=r/1000;
 int j=0;
 std::ofstream out;
 out.open("plot.gp");

 if(save){
   out<<"set terminal pngcairo"<<std::endl;
   out<<"set output 'fig.png'"<<std::endl;
 }
 
 if(cb.name=="EARTH")
   {
     out<<"r = "<<cb.radius<<std::endl;
     out<<"R = "<<cb.radius+1<<std::endl;
     out<<"set title '"<<title<<"'"<<std::endl;
     out<<"#color definitions"<<std::endl;
     out<<"set border lw 1.5"<<std::endl;
     out<<"set style line 1 lc rgb '#000000' lt 1 lw 0.5"<<std::endl;
     out<<"set style line 2 lc rgb '#c0c0c0' lt 2 lw 0.5"<<std::endl;
     out<<"unset key; unset border; set tics scale 0"<<std::endl;
     out<<"set format ''"<<std::endl;
     out<<"set angles degrees"<<std::endl;
     out<<"set xyplane at -1"<<std::endl;
     out<<"set view 75,20"<<std::endl;
     out<<"set lmargin screen 0"<<std::endl;
     out<<"set bmargin screen 0"<<std::endl;
     out<<"set rmargin screen 1"<<std::endl;
     out<<"set tmargin screen 1"<<std::endl;
     out<<"set parametric"<<std::endl;
     out<<"set urange[0:360]"<<std::endl;
     out<<"set vrange[-90:90]"<<std::endl;
     out<<"set isosamples 25"<<std::endl;
     out<<"set hidden3d"<<std::endl;
     out<<"#since we are using Cartesian coordinates, we don't want this"<<std::endl;
     out<<"#set mapping spherical"<<std::endl;
     out<<"splot \
  r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 1,	\
  'world_110m.txt' u (R*cos($1)*cos($2)):(R*sin($1)*cos($2)):(R*sin($2)) w l lw 2 lc rgb 'black'";
 for(auto i: OP){
   out <<", '"<< i.file<<".dat' u 1:2:3 w l lw 3 lc "<<5+j*5;
   j+=1;
 }
 out<<std::endl;
 out<<"pause mouse"<<std::endl;
   }
 else
   {
     out<<"r = "<<cb.radius<<std::endl;
     out<<"R = "<<cb.radius+1<<std::endl;
     out<<"set title '"<<title<<"'"<<std::endl;
     out<<"#color definitions"<<std::endl;
     out<<"set border lw 1.5"<<std::endl;
     out<<"set style line 1 lc rgb '#000000' lt 1 lw 0.5"<<std::endl;
     out<<"set style line 2 lc rgb '#c0c0c0' lt 2 lw 0.5"<<std::endl;
     out<<"unset key; unset border; set tics scale 0"<<std::endl;
     out<<"set format ''"<<std::endl;
     out<<"set angles degrees"<<std::endl;
     out<<"set xyplane at -1"<<std::endl;
     out<<"set view 56,81"<<std::endl;
     out<<"set lmargin screen 0"<<std::endl;
     out<<"set bmargin screen 0"<<std::endl;
     out<<"set rmargin screen 1"<<std::endl;
     out<<"set tmargin screen 1"<<std::endl;
     out<<"set parametric"<<std::endl;
     out<<"set urange[0:360]"<<std::endl;
     out<<"set vrange[-90:90]"<<std::endl;
     out<<"set isosamples 25"<<std::endl;
     out<<"set hidden3d"<<std::endl;
     out<<"#since we are using Cartesian coordinates, we don't want this"<<std::endl;
     out<<"#set mapping spherical"<<std::endl;
     out<<"splot \
  r*cos(v)*cos(u),r*cos(v)*sin(u),r*sin(v) with lines linestyle 1";
     for(auto i: OP){
       out <<", '"<< i.file<<".dat' u 1:2:3 w l lw 3 lc "<<5+j*5;
       j+=1;
     }
     out<<std::endl;
 out<<"pause mouse"<<std::endl;
 
   }
 
 out.close();
 system("gnuplot plot.gp");
}


template <typename T>
std::vector<double> coes2rv(std::vector<double> &coes, bool deg, const T &cb)
{
  double a=0, e=0, i=0, ta=0, aop=0, raan=0, r_norm=0, E=0, rxp=0,ryp=0,rzp=0,vxp=0,vyp=0,vzp=0;
  std::vector <double> rv(6,0);
  std::vector <std::vector <double>> row(3,std::vector<double>(3,0));
  a=coes[0]; e=coes[1]; i=coes[2]; ta=coes[3]; aop=coes[4]; raan=coes[5];
    
  if(deg) //lo pasa a radianes
    {
      i*=d2r;
      ta*=d2r;
      aop*=d2r;
      raan*=d2r;
    }
  
  E=ecc_anomaly(ta,e,"tae"); //calcula la anomalía con newton rapshon
  
  r_norm=a*(1-e*e)/(1+e*std::cos(ta));

  //calcular r y v en el sistema perifocal
  rxp=r_norm*std::cos(ta);  vxp=r_norm*std::sin(E)*std::sqrt(cb.mu*a)/r_norm;
  ryp=r_norm*std::sin(ta);  vyp=std::cos(E)*std::sqrt(1-e*e)*std::sqrt(cb.mu*a)/r_norm;
  rzp=0;                    vzp=0;
 
  std::vector<double> rperif{rxp,ryp,rzp};
  std::vector<double> vperif{vxp,vyp,vzp};
 
  //matriz de rotación 
  row[0]={-std::sin(raan)*std::cos(i)*std::sin(aop)+std::cos(raan)*std::cos(aop),std::cos(raan)*std::cos(i)*std::sin(aop)+std::sin(raan)*std::cos(aop),
	  std::sin(i)*std::sin(aop)};
  row[1]={-std::sin(raan)*std::cos(i)*std::cos(aop)-std::cos(raan)*std::sin(aop),std::cos(raan)*std::cos(i)*std::cos(aop)-std::sin(raan)*std::sin(aop),
	  std::sin(i)*std::cos(aop)};
  row[2]={std::sin(raan)*std::sin(i),-std::cos(raan)*std::sin(i),std::cos(i)};


  //multiclicación de la matriz de rotacion TRANSPUESTA por rperif y vperif
  //para calcular r y v en el sistema inercial del cuerpo central
    
    for(int ii=0; ii<3; ii++)
      for(int jj=0; jj<3; jj++)
	rv[ii]+=rperif[jj]*row[jj][ii];

    
    for(int ii=0; ii<3; ii++)
      for(int jj=0; jj<3; jj++)
	rv[ii+3]+=vperif[jj]*row[jj][ii];

  return rv;
}
double ecc_anomaly(double ta, double e, std::string method)
{
  double  Me=ta, E0=0, E1=0, ratio=0, tol=1e-8;
  
  //newton's method for iteratively finding E
  if (method=="newton")
    {      
      if (Me<M_PI/2.0){
	E0=Me+e/2.0;}
      else{
	E0=Me-e;}
      
      for(int n=0; n<200; n++){
	
	ratio=(E0-e*std::sin(E0)-Me)/(1-e*std::cos(E0));
	if (abs(ratio)<tol)
	  {
	    if (n==0){
	      return E0;}
	    else{
	      return E1;}
	  }
	else{
	  E1=E0-ratio;
	  E0=E1;
	    }
	return 0;
      }
    }
  
  if (method=="tae")
    return 2*std::atan(std::sqrt((1-e)/(1+e))*std::tan(ta/2.0));
  
  else
    std::cout <<"método invalido para calcular la anomalía" << std::endl;
  
  //did not converge
  std::cout << "did not converge" << std::endl;
  return 0;
}

double true_anomaly(double E,double e)
{
return 2*std::atan(std::sqrt((1+e)/(1-e))*std::tan(E/2.0));
}

template <typename T>
std::vector<double> tlecoes(std::string file, const T &cb)
{
  double a,e,i,ta,aop,raan,P,Me,F,E;
  std::string line1, line2, line3;
  std::ifstream data;
  data.open(file);

  std::getline(data,line1);
  std::getline(data,line2);
  std::getline(data,line3);
  
  std::vector<std::string>Line2;
  std::istringstream iss(line2); 
  for(std::string s; iss >> s; ) 
    Line2.push_back(s); 
  
  std::vector<std::string>Line3;
  std::istringstream jss(line3); 
  for(std::string s; jss >> s; ) 
    Line3.push_back(s);
  
  stringstream geeka(Line3[2]); //i
  geeka >> i;
  stringstream geekb(Line3[3]); //raan 
  geekb >> raan; 
  stringstream geekc("0."+Line3[4]); //e
  geekc >> e; 
  stringstream geekd(Line3[5]); //aop 
  geekd >> aop; 
  stringstream geeke(Line3[6]); //Anomalía media 
  geeke >> Me; 
  stringstream geekf(Line3[7]); //frecuencia revs/dia
  geekf >> F;
  P=1/F*24*3600; //segundos
  
  a=std::cbrt((P*P*cb.mu)/(4*std::pow(M_PI,2))); //a

  E=ecc_anomaly(Me, e, "newton");
  
  ta=true_anomaly(E, e);

  std::vector <double> coes {a,e,i,ta,aop,raan};
  return coes;
}
template <typename T>
std::vector <double> rv2coes(std::vector <double> state,const T &CB,bool deg)
{
  double mu=CB.mu; double et=0; double s [6]; double elts[SPICE_OSCLTX_NELTS];

  for(int i=0; i<6; i++)
    s[i]=state[i];

  oscltx_c(s,et,mu,elts); //SPICE FUNCTION:calculate orbital elements for given state
  //these is the output rp,e,i,raan,aop,ma,t0,mu,ta,a,T; in rads
  
  std::vector <double> coes{elts[9],elts[1],elts[2],elts[8],elts[4],elts[3]};
  if(deg)
    {
      coes[2]*=r2d;
      coes[3]*=r2d;
      coes[4]*=r2d;
      coes[5]*=r2d;
    }
  return coes;
}
template <typename T>
void Plot_orbit(std::vector<OrbitPropagator>&OP,T &cb,bool show,bool save,std:: string title,std::vector<std::string> &labels)
  {
    int n=0;
    std::ofstream out;
    std::string save_plot;
    std::string show_plot;

    if(save)
      save_plot="True";
    else
      save_plot="False";
    if(show)
      show_plot="True";
    else
      show_plot="False";
    
    
    out.open("plot.py");
    out<<"import numpy as np"<<std::endl;
    out<<"import plotingfunctions as t"<<std::endl;
    out<<"import planetary_data as pd"<<std::endl;
    out<<"cb=pd."<<cb.name<<std::endl;
    out<<"if __name__ == '__main__':"<<std::endl;
    for(auto i: OP){
      out<<"\tdata"<<n<<"=np.loadtxt('"<<i.file<<".dat',delimiter='\t')"<<std::endl;
      out<<"\tdata"<<n<<"=data"<<n<<"[:,:3]"<<std::endl;
      n+=1;
    }
    n=0;
    out<<"t.plot_n_orbits([";
    
    for(auto i: OP){
      out<<"data"<<n<<",";
      n+=1;
    }
    out<<"],labels=[";
    
    for(auto i: labels)
      out<<"'"<<i<<"',";
    
    out<<"],cb=pd."<<cb.name<<",show_plot="<<show_plot<<",save_plot="<<save_plot<<","<<"title='"<<title<<"')"<<std::endl;
    
    out.close();
    system("python plot.py");
  }
template <typename T>
double calc_atm_density(double z, T &CB)
{
  std::vector <double> data(4,0);
  data=find_rho_z(z,CB);
  if (data[2]==0)
    return 0;
  double Hi=-(data[1]-data[0])/std::log(data[3]/data[2]);
  
  return data[2]*exp(-(z-data[0])/Hi);
}
template <typename T>
std::vector <double> find_rho_z(double z, T &CB)
{
  std::vector <double> f{CB.zs[0],CB.zs[1],CB.zs[2],CB.rhos[0],CB.rhos[1],CB.rhos[2]};
  std::vector <double> ret(4,0);
  if (1>z || z>1000){
    return ret;}
  
  //find the two point surrounding the given input altitude
  for(int n=0; n<2; n++){
    if (f[n]<z<f[n+1])
      {
	ret[0]=f[n]; ret[1]=f[n+1]; ret[2]=f[n+3]; ret[3]=f[n+4];
	return ret;
      }
  }
  return ret;
}
template <typename T>
void Plot_coes(std::vector<OrbitPropagator> &OP,T &CB,bool show, bool save,std::vector<std::string>labels, bool hours, bool days)
{
  calculate_coes(OP,CB); //crea los archivos de coes
  
  std::string save_plot; std::string show_plot; std::string Hours; std::string Days;//parametros de la grafica
  
  int n=0;
  std::ofstream out;
  
  if(save)
    save_plot="True";
  else
    save_plot="False";
  if(show)
    show_plot="True";
  else
    show_plot="False";
  if(hours)
    Hours="True";
  else
    Hours="False";
  if(days)
    Days="True";
  else
    Days="False";
  
  for(auto i: OP){
    std::string s0 = std::to_string(n);
    out.open(s0+"plotcoes.py");
    out<<"import numpy as np"<<std::endl;
    out<<"import plotingfunctions as t"<<std::endl;
    out<<"if __name__ == '__main__':"<<std::endl;
    out<<"\tdata"<<n<<"=np.loadtxt('"<<i.file<<"coes.dat',delimiter='\t')"<<std::endl;
    out<<"t.plot_coes(data"<<n<<",hours="<<Hours<<",days="<<Days<<",show_plot="<<show_plot<<",save_plot="<<
      save_plot<<",title='"<<labels[n]<<" COEs')"<<std::endl;
    n+=1;
    out.close();
  }
  n=0;
  for(auto i: OP)
    {
      std::string s0 = std::to_string(n);
      std::string s1 = "python ";
      std::string s2="plotcoes.py";
      std::string S=s1+s0+s2;
      system((S).c_str());
      n+=1;
    }
}


template <typename T>
void  calculate_coes(std::vector <OrbitPropagator> &OP,T &CB)
{
  std::vector <double> coes (0,6);
  std::vector <double> rv (7,0);
  std::string line;

  for(auto i: OP)
    {
      std::ofstream fout;   //Salida de los datos.                                                                                                                   
      fout.open(i.file+"coes.dat");
      
      std::ifstream data;  //archivo que va a leer
      data.open(i.file+".dat");
      
      while(std::getline(data,line))
	{
	  std::vector<std::string>Line;
	  std::istringstream iss(line); 
	  for(std::string s; iss >> s; ) 
	    Line.push_back(s);
      
	  for(int i=0; i<7; i++) //copiar la linea de rv en un vector
	    {
	      stringstream geeka(Line[i]); 
	      geeka >> rv[i];
	    }
	  coes=rv2coes(rv,CB,true); //el tamaño de rv no importa mientras los 6 primeros datos sean x,y,z,vx,vy,vz
  
	  for(auto j: coes)
	    fout <<j<<"\t"; //estos son los coes
	  fout<<rv[6]<<std::endl; //este es el tiempo
	}
       fout.close();
       data.close();
    }
}
