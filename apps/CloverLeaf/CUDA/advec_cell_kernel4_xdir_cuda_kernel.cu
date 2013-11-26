//
// auto-generated by ops.py on 2013-11-26 11:43
//

__constant__ int xdim0_advec_cell_kernel4_xdir;
__constant__ int xdim1_advec_cell_kernel4_xdir;
__constant__ int xdim2_advec_cell_kernel4_xdir;
__constant__ int xdim3_advec_cell_kernel4_xdir;
__constant__ int xdim4_advec_cell_kernel4_xdir;
__constant__ int xdim5_advec_cell_kernel4_xdir;
__constant__ int xdim6_advec_cell_kernel4_xdir;
__constant__ int xdim7_advec_cell_kernel4_xdir;
__constant__ int xdim8_advec_cell_kernel4_xdir;
__constant__ int xdim9_advec_cell_kernel4_xdir;
__constant__ int xdim10_advec_cell_kernel4_xdir;

#define OPS_ACC0(x,y) (x+xdim0_advec_cell_kernel4_xdir*(y))
#define OPS_ACC1(x,y) (x+xdim1_advec_cell_kernel4_xdir*(y))
#define OPS_ACC2(x,y) (x+xdim2_advec_cell_kernel4_xdir*(y))
#define OPS_ACC3(x,y) (x+xdim3_advec_cell_kernel4_xdir*(y))
#define OPS_ACC4(x,y) (x+xdim4_advec_cell_kernel4_xdir*(y))
#define OPS_ACC5(x,y) (x+xdim5_advec_cell_kernel4_xdir*(y))
#define OPS_ACC6(x,y) (x+xdim6_advec_cell_kernel4_xdir*(y))
#define OPS_ACC7(x,y) (x+xdim7_advec_cell_kernel4_xdir*(y))
#define OPS_ACC8(x,y) (x+xdim8_advec_cell_kernel4_xdir*(y))
#define OPS_ACC9(x,y) (x+xdim9_advec_cell_kernel4_xdir*(y))
#define OPS_ACC10(x,y) (x+xdim10_advec_cell_kernel4_xdir*(y))

//user function
__device__

inline void advec_cell_kernel4_xdir( double *density1, double *energy1,
                         const double *mass_flux_x, const double *vol_flux_x,
                         double *pre_vol, double *post_vol,
                         double *pre_mass, double *post_mass,
                         double *advec_vol, double *post_ener,
                         double *ener_flux) {

  pre_mass[OPS_ACC6(0,0)] = density1[OPS_ACC0(0,0)] * pre_vol[OPS_ACC4(0,0)];
  post_mass[OPS_ACC7(0,0)] = pre_mass[OPS_ACC6(0,0)] + mass_flux_x[OPS_ACC2(0,0)] - mass_flux_x[OPS_ACC2(1,0)];
  post_ener[OPS_ACC9(0,0)] = ( energy1[OPS_ACC1(0,0)] * pre_mass[OPS_ACC6(0,0)] + ener_flux[OPS_ACC10(0,0)] - ener_flux[OPS_ACC10(1,0)])/post_mass[OPS_ACC7(0,0)];
  advec_vol[OPS_ACC8(0,0)] = pre_vol[OPS_ACC4(0,0)] + vol_flux_x[OPS_ACC3(0,0)] - vol_flux_x[OPS_ACC3(1,0)];
  density1[OPS_ACC0(0,0)] = post_mass[OPS_ACC7(0,0)]/advec_vol[OPS_ACC8(0,0)];
  energy1[OPS_ACC1(0,0)] = post_ener[OPS_ACC9(0,0)];

}



#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6
#undef OPS_ACC7
#undef OPS_ACC8
#undef OPS_ACC9
#undef OPS_ACC10


__global__ void ops_advec_cell_kernel4_xdir(
double* __restrict arg0,
double* __restrict arg1,
const double* __restrict arg2,
const double* __restrict arg3,
double* __restrict arg4,
double* __restrict arg5,
double* __restrict arg6,
double* __restrict arg7,
double* __restrict arg8,
double* __restrict arg9,
double* __restrict arg10,
int size0,
int size1 ){


  int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1 + idx_y * 1 * xdim0_advec_cell_kernel4_xdir;
  arg1 += idx_x * 1 + idx_y * 1 * xdim1_advec_cell_kernel4_xdir;
  arg2 += idx_x * 1 + idx_y * 1 * xdim2_advec_cell_kernel4_xdir;
  arg3 += idx_x * 1 + idx_y * 1 * xdim3_advec_cell_kernel4_xdir;
  arg4 += idx_x * 1 + idx_y * 1 * xdim4_advec_cell_kernel4_xdir;
  arg5 += idx_x * 1 + idx_y * 1 * xdim5_advec_cell_kernel4_xdir;
  arg6 += idx_x * 1 + idx_y * 1 * xdim6_advec_cell_kernel4_xdir;
  arg7 += idx_x * 1 + idx_y * 1 * xdim7_advec_cell_kernel4_xdir;
  arg8 += idx_x * 1 + idx_y * 1 * xdim8_advec_cell_kernel4_xdir;
  arg9 += idx_x * 1 + idx_y * 1 * xdim9_advec_cell_kernel4_xdir;
  arg10 += idx_x * 1 + idx_y * 1 * xdim10_advec_cell_kernel4_xdir;

  if (idx_x < size0 && idx_y < size1) {
    advec_cell_kernel4_xdir(arg0, arg1, arg2, arg3,
                   arg4, arg5, arg6, arg7, arg8,
                   arg9, arg10);
  }

}

// host stub function
void ops_par_loop_advec_cell_kernel4_xdir(char const *name, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6, ops_arg arg7, ops_arg arg8,
 ops_arg arg9, ops_arg arg10) {

  ops_arg args[11] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10};


  int x_size = range[1]-range[0];
  int y_size = range[3]-range[2];

  int xdim0 = args[0].dat->block_size[0];
  int xdim1 = args[1].dat->block_size[0];
  int xdim2 = args[2].dat->block_size[0];
  int xdim3 = args[3].dat->block_size[0];
  int xdim4 = args[4].dat->block_size[0];
  int xdim5 = args[5].dat->block_size[0];
  int xdim6 = args[6].dat->block_size[0];
  int xdim7 = args[7].dat->block_size[0];
  int xdim8 = args[8].dat->block_size[0];
  int xdim9 = args[9].dat->block_size[0];
  int xdim10 = args[10].dat->block_size[0];


  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(8,"advec_cell_kernel4_xdir");
  ops_timers_core(&c1,&t1);

  if (OPS_kernels[8].count == 0) {
    cudaMemcpyToSymbol( xdim0_advec_cell_kernel4_xdir, &xdim0, sizeof(int) );
    cudaMemcpyToSymbol( xdim1_advec_cell_kernel4_xdir, &xdim1, sizeof(int) );
    cudaMemcpyToSymbol( xdim2_advec_cell_kernel4_xdir, &xdim2, sizeof(int) );
    cudaMemcpyToSymbol( xdim3_advec_cell_kernel4_xdir, &xdim3, sizeof(int) );
    cudaMemcpyToSymbol( xdim4_advec_cell_kernel4_xdir, &xdim4, sizeof(int) );
    cudaMemcpyToSymbol( xdim5_advec_cell_kernel4_xdir, &xdim5, sizeof(int) );
    cudaMemcpyToSymbol( xdim6_advec_cell_kernel4_xdir, &xdim6, sizeof(int) );
    cudaMemcpyToSymbol( xdim7_advec_cell_kernel4_xdir, &xdim7, sizeof(int) );
    cudaMemcpyToSymbol( xdim8_advec_cell_kernel4_xdir, &xdim8, sizeof(int) );
    cudaMemcpyToSymbol( xdim9_advec_cell_kernel4_xdir, &xdim9, sizeof(int) );
    cudaMemcpyToSymbol( xdim10_advec_cell_kernel4_xdir, &xdim10, sizeof(int) );
  }



  dim3 grid( (x_size-1)/OPS_block_size_x+ 1, (y_size-1)/OPS_block_size_y + 1, 1);
  dim3 block(OPS_block_size_x,OPS_block_size_y,1);




  char *p_a[11];


  //set up initial pointers
  p_a[0] = &args[0].data_d[
  + args[0].dat->size * args[0].dat->block_size[0] * ( range[2] * 1 - args[0].dat->offset[1] )
  + args[0].dat->size * ( range[0] * 1 - args[0].dat->offset[0] ) ];

  p_a[1] = &args[1].data_d[
  + args[1].dat->size * args[1].dat->block_size[0] * ( range[2] * 1 - args[1].dat->offset[1] )
  + args[1].dat->size * ( range[0] * 1 - args[1].dat->offset[0] ) ];

  p_a[2] = &args[2].data_d[
  + args[2].dat->size * args[2].dat->block_size[0] * ( range[2] * 1 - args[2].dat->offset[1] )
  + args[2].dat->size * ( range[0] * 1 - args[2].dat->offset[0] ) ];

  p_a[3] = &args[3].data_d[
  + args[3].dat->size * args[3].dat->block_size[0] * ( range[2] * 1 - args[3].dat->offset[1] )
  + args[3].dat->size * ( range[0] * 1 - args[3].dat->offset[0] ) ];

  p_a[4] = &args[4].data_d[
  + args[4].dat->size * args[4].dat->block_size[0] * ( range[2] * 1 - args[4].dat->offset[1] )
  + args[4].dat->size * ( range[0] * 1 - args[4].dat->offset[0] ) ];

  p_a[5] = &args[5].data_d[
  + args[5].dat->size * args[5].dat->block_size[0] * ( range[2] * 1 - args[5].dat->offset[1] )
  + args[5].dat->size * ( range[0] * 1 - args[5].dat->offset[0] ) ];

  p_a[6] = &args[6].data_d[
  + args[6].dat->size * args[6].dat->block_size[0] * ( range[2] * 1 - args[6].dat->offset[1] )
  + args[6].dat->size * ( range[0] * 1 - args[6].dat->offset[0] ) ];

  p_a[7] = &args[7].data_d[
  + args[7].dat->size * args[7].dat->block_size[0] * ( range[2] * 1 - args[7].dat->offset[1] )
  + args[7].dat->size * ( range[0] * 1 - args[7].dat->offset[0] ) ];

  p_a[8] = &args[8].data_d[
  + args[8].dat->size * args[8].dat->block_size[0] * ( range[2] * 1 - args[8].dat->offset[1] )
  + args[8].dat->size * ( range[0] * 1 - args[8].dat->offset[0] ) ];

  p_a[9] = &args[9].data_d[
  + args[9].dat->size * args[9].dat->block_size[0] * ( range[2] * 1 - args[9].dat->offset[1] )
  + args[9].dat->size * ( range[0] * 1 - args[9].dat->offset[0] ) ];

  p_a[10] = &args[10].data_d[
  + args[10].dat->size * args[10].dat->block_size[0] * ( range[2] * 1 - args[10].dat->offset[1] )
  + args[10].dat->size * ( range[0] * 1 - args[10].dat->offset[0] ) ];


  ops_halo_exchanges_cuda(args, 11);


  //call kernel wrapper function, passing in pointers to data
  ops_advec_cell_kernel4_xdir<<<grid, block >>> (  (double *)p_a[0], (double *)p_a[1],
           (double *)p_a[2], (double *)p_a[3],
           (double *)p_a[4], (double *)p_a[5],
           (double *)p_a[6], (double *)p_a[7],
           (double *)p_a[8], (double *)p_a[9],
           (double *)p_a[10],x_size, y_size);

  if (OPS_diags>1) cudaDeviceSynchronize();
  ops_set_dirtybit_cuda(args, 11);

  //Update kernel record
  ops_timers_core(&c2,&t2);
  OPS_kernels[8].count++;
  OPS_kernels[8].time += t2-t1;
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg2);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg3);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg4);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg5);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg6);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg7);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg8);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg9);
  OPS_kernels[8].transfer += ops_compute_transfer(dim, range, &arg10);
}
