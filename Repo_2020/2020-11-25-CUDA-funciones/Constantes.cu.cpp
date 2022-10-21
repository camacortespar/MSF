//Mi Primer Programa en CUDA
//scp MiPrimerProgramaEnCuda.cu cdelv@168.176.8.34:CUDA2020
//nvcc -arch=sm_30  para darle un nombre -o ArchivoCUDA
//time ./ExplosiveSOC_LIST_SSF

#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Q 5
/*
#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;
*/
//--------------- KERNELS ----------------
__constant__ float d_w[Q];
__constant__ int d_Vy[Q];
__constant__ int d_Vx[Q];

__global__ void Test(float *d_test){
  d_test[0]=d_Vx[0];
}

int main(){
  //DECLARAR LAS MATRICES
  float h_w[Q];
  int h_Vy[Q];
  int h_Vx[Q];
  float h_test[1]; float*d_test; cudaMalloc((void**) &d_test,sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  h_w[0]=1.0/3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  h_Vx[0]=0, h_Vx[1]=1, h_Vx[2]=0,  h_Vx[3]=-1,  h_Vx[4]=0;
  h_Vy[0]=0, h_Vy[1]=0, h_Vy[2]=1,  h_Vy[3]=0,   h_Vy[4]=-1;
 

  
  //Enviarlos al Device
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(float),0,cudaMemcpyHostToDevice);

  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(1,1,1);
  dim3 BlocksPerGrid(1,1,1);
  Test<<<BlocksPerGrid,ThreadsPerBlock>>>(d_test);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_test,d_test,sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  cout<<h_test[0]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_test);

  return 0;
}
