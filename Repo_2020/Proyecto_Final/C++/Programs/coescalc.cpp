#include "master.h"


int main(void)
{
  earth cb;
  
  std::vector<double> coes{6,0};
  std::vector<double> rv{6791.72590054471, 230.74102783728756, 1.9488697305752931, -0.2600812919334669, 7.653674225886742, 0.12046463959766343};

  coes=rv2coes(rv,cb,true);

  std::cout <<"a="<<coes[0]<<std::endl;
  std::cout <<"e="<<coes[1]<<std::endl;
  std::cout <<"i="<<coes[2]<<std::endl;
  std::cout <<"ta="<<coes[3]<<std::endl;
  std::cout <<"aop="<<coes[4]<<std::endl;
  std::cout <<"raan="<<coes[5]<<std::endl;
  
  return 0;
}
