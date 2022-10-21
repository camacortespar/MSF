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
#define Ly 8
#define Nx 8
#define Ny 8
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Ny;

//--------------- KERNELS ----------------
__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];

__device__ float SumeleUno(float x){
  return x+1;
}
__global__ void SumarleUnoATodos(float * d_a,size_t pitcha){
  int ix,iy; float *aux;
  ix=blockIdx.x*blockDim.x+threadIdx.x;
  iy=blockIdx.y*blockDim.y+threadIdx.y;

  aux=d_a+(ix*pitcha)/sizeof(float)+iy; // aux es &(d_a[ix][iy])
  
  (*aux)=SumeleUno(*aux); //  (*aux) es d_a[ix][iy]
}

int main(){
  int ix,iy;
  //DECLARAR LAS MATRICES
  float h_w[Q]; int h_Vx[Q],h_Vy[Q];

  //INICIALIZAR LAS CONSTANTES Y LAS MANDO AL DEVICE
  //Cargarlos en el Host
  h_w[0]=1.0/3; h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  //Enviarlos al Device
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice);

  //DECLARAR LAS VARIABLES
  //declarar en Host
  float h_a[Lx][Ly];
  //declarar en Device
  float*d_a; size_t pitcha;
  cudaMallocPitch((void**) &d_a,&pitcha,Ly*sizeof(float),Lx);

  //INICIALIZAR LAS VARIABLES Y LAS MANDO AL DEVICE
  //Cargar los datos en el Host
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      h_a[ix][iy]=Ly*ix+iy;
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;
  //Enviarlos al Device
  cudaMemcpy2D(d_a,pitcha,h_a,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);

  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  SumarleUnoATodos<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,pitcha);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy2D(h_a,Ly*sizeof(float),d_a,pitcha,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
