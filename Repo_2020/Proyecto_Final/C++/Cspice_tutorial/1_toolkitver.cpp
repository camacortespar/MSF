#include <iostream>
#include <string>
#include <stdio.h>
#include <SpiceUsr.h>
 
int main()
{
  std::string versn;
 
  versn = tkvrsn_c( "TOOLKIT" );
 
  std::cout<< "Toolkit version "<< versn <<std::endl;;
 
  return(0);
}
