//
// auto-generated by ops.py
//
#include "./OpenACC/poisson_common.h"

#define OPS_GPU

int xdim3_poisson_kernel_populate;
int xdim4_poisson_kernel_populate;
int xdim5_poisson_kernel_populate;

#define OPS_ACC3(x,y) (x+xdim3_poisson_kernel_populate*(y))
#define OPS_ACC4(x,y) (x+xdim4_poisson_kernel_populate*(y))
#define OPS_ACC5(x,y) (x+xdim5_poisson_kernel_populate*(y))

//user function
inline
void poisson_kernel_populate(const int *dispx, const int *dispy, const int *idx, double *u, double *f, double *ref) {
  double x = dx * (double)(idx[0]+dispx[0]);
  double y = dy * (double)(idx[1]+dispy[0]);

  u[OPS_ACC3(0,0)] = sin(M_PI*x)*cos(2.0*M_PI*y);
  f[OPS_ACC4(0,0)] = -5.0*M_PI*M_PI*sin(M_PI*x)*cos(2.0*M_PI*y);
  ref[OPS_ACC5(0,0)] = sin(M_PI*x)*cos(2.0*M_PI*y);

}


#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5



void poisson_kernel_populate_c_wrapper(
  int p_a0,
  int p_a1,
  int *p_a2,
  double *p_a3,
  double *p_a4,
  double *p_a5,
  int arg_idx0, int arg_idx1,
  int x_size, int y_size) {
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a3,p_a4,p_a5)
  #pragma acc loop
  #endif
  for ( int n_y=0; n_y<y_size; n_y++ ){
    #ifdef OPS_GPU
    #pragma acc loop
    #endif
    for ( int n_x=0; n_x<x_size; n_x++ ){
      int arg_idx[] = {arg_idx0+n_x, arg_idx1+n_y};
      poisson_kernel_populate(  &p_a0,
           &p_a1,arg_idx,
           p_a3 + n_x*1*1 + n_y*xdim3_poisson_kernel_populate*1*1, p_a4 + n_x*1*1 + n_y*xdim4_poisson_kernel_populate*1*1,
           p_a5 + n_x*1*1 + n_y*xdim5_poisson_kernel_populate*1*1 );

    }
  }
}
