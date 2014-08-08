//
// auto-generated by ops.py
//

#ifdef OCL_FMA
#pragma OPENCL FP_CONTRACT ON
#else
#pragma OPENCL FP_CONTRACT OFF
#endif
#pragma OPENCL EXTENSION cl_khr_fp64:enable

#include "user_types.h"
#include "ops_opencl_reduction.h"

#ifndef MIN
#define MIN(a,b) ((a<b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a>b) ? (a) : (b))
#endif
#ifndef SIGN
#define SIGN(a,b) ((b<0.0) ? (a*(-1)) : (a))
#endif
#define OPS_READ 0
#define OPS_WRITE 1
#define OPS_RW 2
#define OPS_INC 3
#define OPS_MIN 4
#define OPS_MAX 5
#define ZERO_double 0.0;
#define INFINITY_double INFINITY;
#define ZERO_float 0.0f;
#define INFINITY_float INFINITY;
#define ZERO_int 0;
#define INFINITY_int INFINITY;
#define ZERO_uint 0;
#define INFINITY_uint INFINITY;
#define ZERO_ll 0;
#define INFINITY_ll INFINITY;
#define ZERO_ull 0;
#define INFINITY_ull INFINITY;
#define ZERO_bool 0;
#define OPS_ACC0(x,y,z) (x+xdim0_reset_field_kernel2*(y)+xdim0_reset_field_kernel2*ydim0_reset_field_kernel2*(z))
#define OPS_ACC1(x,y,z) (x+xdim1_reset_field_kernel2*(y)+xdim1_reset_field_kernel2*ydim1_reset_field_kernel2*(z))
#define OPS_ACC2(x,y,z) (x+xdim2_reset_field_kernel2*(y)+xdim2_reset_field_kernel2*ydim2_reset_field_kernel2*(z))
#define OPS_ACC3(x,y,z) (x+xdim3_reset_field_kernel2*(y)+xdim3_reset_field_kernel2*ydim3_reset_field_kernel2*(z))
#define OPS_ACC4(x,y,z) (x+xdim4_reset_field_kernel2*(y)+xdim4_reset_field_kernel2*ydim4_reset_field_kernel2*(z))
#define OPS_ACC5(x,y,z) (x+xdim5_reset_field_kernel2*(y)+xdim5_reset_field_kernel2*ydim5_reset_field_kernel2*(z))


//user function
void reset_field_kernel2( __global double * restrict xvel0, const __global double * restrict xvel1, __global double * restrict yvel0, 
const __global double * restrict yvel1, __global double * restrict zvel0, const __global double * restrict zvel1)

  {

  xvel0[OPS_ACC0(0,0,0)]  = xvel1[OPS_ACC1(0,0,0)] ;
  yvel0[OPS_ACC2(0,0,0)]  = yvel1[OPS_ACC3(0,0,0)] ;
  zvel0[OPS_ACC4(0,0,0)]  = zvel1[OPS_ACC5(0,0,0)] ;
}



 #undef OPS_ACC0
 #undef OPS_ACC1
 #undef OPS_ACC2
 #undef OPS_ACC3
 #undef OPS_ACC4
 #undef OPS_ACC5


 __kernel void ops_reset_field_kernel2(
 __global double* restrict arg0,
 __global const double* restrict arg1,
 __global double* restrict arg2,
 __global const double* restrict arg3,
 __global double* restrict arg4,
 __global const double* restrict arg5,
 const int base0,
 const int base1,
 const int base2,
 const int base3,
 const int base4,
 const int base5,
 const int size0,
 const int size1,
 const int size2 ){


   int idx_z = get_global_id(2);
   int idx_y = get_global_id(1);
   int idx_x = get_global_id(0);

   if (idx_x < size0 && idx_y < size1 && idx_z < size2) {
     reset_field_kernel2(&arg0[base0 + idx_x * 1 + idx_y * 1 * xdim0_reset_field_kernel2 + idx_z * 1 * xdim0_reset_field_kernel2 * ydim0_reset_field_kernel2],
                        &arg1[base1 + idx_x * 1 + idx_y * 1 * xdim1_reset_field_kernel2 + idx_z * 1 * xdim1_reset_field_kernel2 * ydim1_reset_field_kernel2],
                        &arg2[base2 + idx_x * 1 + idx_y * 1 * xdim2_reset_field_kernel2 + idx_z * 1 * xdim2_reset_field_kernel2 * ydim2_reset_field_kernel2],
                        &arg3[base3 + idx_x * 1 + idx_y * 1 * xdim3_reset_field_kernel2 + idx_z * 1 * xdim3_reset_field_kernel2 * ydim3_reset_field_kernel2],
                        &arg4[base4 + idx_x * 1 + idx_y * 1 * xdim4_reset_field_kernel2 + idx_z * 1 * xdim4_reset_field_kernel2 * ydim4_reset_field_kernel2],
                        &arg5[base5 + idx_x * 1 + idx_y * 1 * xdim5_reset_field_kernel2 + idx_z * 1 * xdim5_reset_field_kernel2 * ydim5_reset_field_kernel2]);
   }

 }
