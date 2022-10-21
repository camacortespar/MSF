//Mi Primer Programa en CUDA
#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

//--------------- KERNELS ----------------
__global__ void Test(float * d_test){
  (*d_test)=5.0;
}

int main(){
  //DECLARAR LAS MATRICES
  float h_test[1];
  float *d_test;  cudaMalloc((void**) &d_test,sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  h_test[0]=2.5;
  //Enviarlos al Device
  cudaMemcpy(d_test,h_test,sizeof(float),cudaMemcpyHostToDevice);
  
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
