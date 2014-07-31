//
// auto-generated by ops.py on 2014-07-31 12:19
//


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
#define OPS_ACC0(x,y,z) (x+xdim0_field_summary_kernel*(y)+xdim0_field_summary_kernel*ydim0_field_summary_kernel*(z))
#define OPS_ACC1(x,y,z) (x+xdim1_field_summary_kernel*(y)+xdim1_field_summary_kernel*ydim1_field_summary_kernel*(z))
#define OPS_ACC2(x,y,z) (x+xdim2_field_summary_kernel*(y)+xdim2_field_summary_kernel*ydim2_field_summary_kernel*(z))
#define OPS_ACC3(x,y,z) (x+xdim3_field_summary_kernel*(y)+xdim3_field_summary_kernel*ydim3_field_summary_kernel*(z))
#define OPS_ACC4(x,y,z) (x+xdim4_field_summary_kernel*(y)+xdim4_field_summary_kernel*ydim4_field_summary_kernel*(z))
#define OPS_ACC5(x,y,z) (x+xdim5_field_summary_kernel*(y)+xdim5_field_summary_kernel*ydim5_field_summary_kernel*(z))
#define OPS_ACC6(x,y,z) (x+xdim6_field_summary_kernel*(y)+xdim6_field_summary_kernel*ydim6_field_summary_kernel*(z))


//user function
void field_summary_kernel( const __global double * restrict volume, const __global double * restrict density0, const __global double * restrict energy0, 
const __global double * restrict pressure, const __global double * restrict xvel0, const __global double * restrict yvel0, const __global double * restrict zvel0, 
 double * restrict vol, double * restrict mass, double * restrict ie, double * restrict ke,
 double * restrict press)

  {

  double vsqrd, cell_vol, cell_mass;

  vsqrd = 0.0;
  for (int iz=0;iz<2;iz++){
    for (int iy=0;iy<2;iy++){
      for (int ix=0;ix<2;ix++){
        vsqrd+=0.125*( xvel0[OPS_ACC4(ix,iy,iz)] * xvel0[OPS_ACC4(ix,iy,iz)] +
                       yvel0[OPS_ACC5(ix,iy,iz)] * yvel0[OPS_ACC5(ix,iy,iz)] +
                       zvel0[OPS_ACC6(ix,iy,iz)] * zvel0[OPS_ACC6(ix,iy,iz)]);
      }
    }
  }

  cell_vol = volume[OPS_ACC0(0,0,0)];
  cell_mass = cell_vol * density0[OPS_ACC1(0,0,0)];
  *vol = *vol + cell_vol;
  *mass = *mass + cell_mass;
  *ie = *ie + cell_mass * energy0[OPS_ACC2(0,0,0)];
  *ke = *ke + cell_mass * 0.5 * vsqrd;
  *press = *press + cell_vol * pressure[OPS_ACC3(0,0,0)];

}



 #undef OPS_ACC0
 #undef OPS_ACC1
 #undef OPS_ACC2
 #undef OPS_ACC3
 #undef OPS_ACC4
 #undef OPS_ACC5
 #undef OPS_ACC6


 __kernel void ops_field_summary_kernel(
 __global const double* restrict arg0,
 __global const double* restrict arg1,
 __global const double* restrict arg2,
 __global const double* restrict arg3,
 __global const double* restrict arg4,
 __global const double* restrict arg5,
 __global const double* restrict arg6,
 __global double* restrict arg7,
 __local double* scratch7,
 int r_bytes7,
 __global double* restrict arg8,
 __local double* scratch8,
 int r_bytes8,
 __global double* restrict arg9,
 __local double* scratch9,
 int r_bytes9,
 __global double* restrict arg10,
 __local double* scratch10,
 int r_bytes10,
 __global double* restrict arg11,
 __local double* scratch11,
 int r_bytes11,
 const int base0,
 const int base1,
 const int base2,
 const int base3,
 const int base4,
 const int base5,
 const int base6,
 const int size0,
 const int size1,
 const int size2 ){

   arg7 += r_bytes7;
   double arg7_l[1];
   arg8 += r_bytes8;
   double arg8_l[1];
   arg9 += r_bytes9;
   double arg9_l[1];
   arg10 += r_bytes10;
   double arg10_l[1];
   arg11 += r_bytes11;
   double arg11_l[1];
   for (int d=0; d<1; d++) arg7_l[d] = ZERO_double;
   for (int d=0; d<1; d++) arg8_l[d] = ZERO_double;
   for (int d=0; d<1; d++) arg9_l[d] = ZERO_double;
   for (int d=0; d<1; d++) arg10_l[d] = ZERO_double;
   for (int d=0; d<1; d++) arg11_l[d] = ZERO_double;

   int idx_z = get_global_id(2);
   int idx_y = get_global_id(1);
   int idx_x = get_global_id(0);

   if (idx_x < size0 && idx_y < size1 && idx_z < size2) {
     field_summary_kernel(&arg0[base0 + idx_x * 1 + idx_y * 1 * xdim0_field_summary_kernel + idx_z * 1 * xdim0_field_summary_kernel * ydim0_field_summary_kernel],
                          &arg1[base1 + idx_x * 1 + idx_y * 1 * xdim1_field_summary_kernel + idx_z * 1 * xdim1_field_summary_kernel * ydim1_field_summary_kernel],
                          &arg2[base2 + idx_x * 1 + idx_y * 1 * xdim2_field_summary_kernel + idx_z * 1 * xdim2_field_summary_kernel * ydim2_field_summary_kernel],
                          &arg3[base3 + idx_x * 1 + idx_y * 1 * xdim3_field_summary_kernel + idx_z * 1 * xdim3_field_summary_kernel * ydim3_field_summary_kernel],
                          &arg4[base4 + idx_x * 1 + idx_y * 1 * xdim4_field_summary_kernel + idx_z * 1 * xdim4_field_summary_kernel * ydim4_field_summary_kernel],
                          &arg5[base5 + idx_x * 1 + idx_y * 1 * xdim5_field_summary_kernel + idx_z * 1 * xdim5_field_summary_kernel * ydim5_field_summary_kernel],
                          &arg6[base6 + idx_x * 1 + idx_y * 1 * xdim6_field_summary_kernel + idx_z * 1 * xdim6_field_summary_kernel * ydim6_field_summary_kernel],
                          arg7_l,
                          arg8_l,
                          arg9_l,
                          arg10_l,
                          arg11_l);
   }
   reduce_double(arg7_l[0], scratch7, arg7, OPS_INC);
   reduce_double(arg8_l[0], scratch8, arg8, OPS_INC);
   reduce_double(arg9_l[0], scratch9, arg9, OPS_INC);
   reduce_double(arg10_l[0], scratch10, arg10, OPS_INC);
   reduce_double(arg11_l[0], scratch11, arg11, OPS_INC);

 }
