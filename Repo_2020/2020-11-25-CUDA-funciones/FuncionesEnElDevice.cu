//Mi Primer Programa en CUDA
#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Q 5

#define Lx 16
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------- KERNELS ----------------
__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];

__device__ float SumeleUno(float x){
  return x+1;
}
__global__ void SumarleUnoATodos(float * d_a){
  int ix;
  ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_a[ix]=SumeleUno(d_a[ix]);
}

int main(){
  int ix;
  //DECLARAR LAS MATRICES
  float h_w[Q]; int h_Vx[Q],h_Vy[Q];
  float h_a[Lx]; float*d_a; cudaMalloc((void**) &d_a,Lx*sizeof(float));

  //INICIALIZAR LAS CONSTANTES Y LAS MANDO AL DEVICE
  //Cargarlos en el Host
  h_w[0]=1.0/3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  //Enviarlos al Device
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice);

  //INICIALIZAR LAS VARIABLES Y LAS MANDO AL DEVICE
  for(ix=0;ix<Lx;ix++) h_a[ix]=ix;
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);

  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumarleUnoATodos<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_a,d_a,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<ix<<" "<<h_a[ix]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
