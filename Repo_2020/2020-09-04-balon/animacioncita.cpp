//entubar esto con gnuplot para hacer una animacion
//./a.out | gnuplot
#include <iostream>
#include <cmath>


 
int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);


  std::cout <<"plot \"datos.dat\" u 2:3 w l" << std::endl;
  std::cout <<"pause 3" << std::endl;
 
  return 0;
}
