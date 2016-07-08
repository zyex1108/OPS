//
// auto-generated by ops.py
//
#include "./OpenACC/clover_leaf_common.h"

#undef OPS_GPU

int xdim0_reset_field_kernel1;
int ydim0_reset_field_kernel1;
int xdim1_reset_field_kernel1;
int ydim1_reset_field_kernel1;
int xdim2_reset_field_kernel1;
int ydim2_reset_field_kernel1;
int xdim3_reset_field_kernel1;
int ydim3_reset_field_kernel1;


#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3


#define OPS_ACC0(x,y,z) (x+xdim0_reset_field_kernel1*(y)+xdim0_reset_field_kernel1*ydim0_reset_field_kernel1*(z))
#define OPS_ACC1(x,y,z) (x+xdim1_reset_field_kernel1*(y)+xdim1_reset_field_kernel1*ydim1_reset_field_kernel1*(z))
#define OPS_ACC2(x,y,z) (x+xdim2_reset_field_kernel1*(y)+xdim2_reset_field_kernel1*ydim2_reset_field_kernel1*(z))
#define OPS_ACC3(x,y,z) (x+xdim3_reset_field_kernel1*(y)+xdim3_reset_field_kernel1*ydim3_reset_field_kernel1*(z))

//user function
inline
void reset_field_kernel1( double *density0, const double *density1,
                          double *energy0,  const double *energy1) {

  density0[OPS_ACC0(0,0,0)]  = density1[OPS_ACC1(0,0,0)] ;
  energy0[OPS_ACC2(0,0,0)]  = energy1[OPS_ACC3(0,0,0)] ;

}


#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3



void reset_field_kernel1_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  double *p_a3,
  int x_size, int y_size, int z_size) {
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a0,p_a1,p_a2,p_a3)
  #pragma acc loop
  #endif
  for ( int n_z=0; n_z<z_size; n_z++ ){
    #ifdef OPS_GPU
    #pragma acc loop
    #endif
    for ( int n_y=0; n_y<y_size; n_y++ ){
      #ifdef OPS_GPU
      #pragma acc loop
      #endif
      for ( int n_x=0; n_x<x_size; n_x++ ){
        reset_field_kernel1(  p_a0 + n_x*1*1 + n_y*xdim0_reset_field_kernel1*1*1 + n_z*xdim0_reset_field_kernel1*ydim0_reset_field_kernel1*1,
           p_a1 + n_x*1*1 + n_y*xdim1_reset_field_kernel1*1*1 + n_z*xdim1_reset_field_kernel1*ydim1_reset_field_kernel1*1,
           p_a2 + n_x*1*1 + n_y*xdim2_reset_field_kernel1*1*1 + n_z*xdim2_reset_field_kernel1*ydim2_reset_field_kernel1*1,
           p_a3 + n_x*1*1 + n_y*xdim3_reset_field_kernel1*1*1 + n_z*xdim3_reset_field_kernel1*ydim3_reset_field_kernel1*1 );

      }
    }
  }
}
