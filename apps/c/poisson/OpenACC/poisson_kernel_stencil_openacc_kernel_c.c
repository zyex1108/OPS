//
// auto-generated by ops.py
//
#include "./OpenACC/poisson_common.h"

#undef OPS_GPU

int xdim0_poisson_kernel_stencil;
int xdim1_poisson_kernel_stencil;


#undef OPS_ACC0
#undef OPS_ACC1


#define OPS_ACC0(x,y) (x+xdim0_poisson_kernel_stencil*(y))
#define OPS_ACC1(x,y) (x+xdim1_poisson_kernel_stencil*(y))

//user function
inline
void poisson_kernel_stencil(const double *u, double *u2) {
  u2[OPS_ACC1(0,0)] = ((u[OPS_ACC0(-1,0)]+u[OPS_ACC0(1,0)])*0.125f
                     + (u[OPS_ACC0(0,-1)]+u[OPS_ACC0(0,1)])*0.125f
                     - 0.125f*u[OPS_ACC0(0,0)]);
}


#undef OPS_ACC0
#undef OPS_ACC1



void poisson_kernel_stencil_c_wrapper(
  double *p_a0,
  double *p_a1,
  int x_size, int y_size) {
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a0,p_a1)
  #pragma acc loop
  #endif
  for ( int n_y=0; n_y<y_size; n_y++ ){
    #ifdef OPS_GPU
    #pragma acc loop
    #endif
    for ( int n_x=0; n_x<x_size; n_x++ ){
      poisson_kernel_stencil(  p_a0 + n_x*1*1 + n_y*xdim0_poisson_kernel_stencil*1*1,
           p_a1 + n_x*1*1 + n_y*xdim1_poisson_kernel_stencil*1*1 );

    }
  }
}
