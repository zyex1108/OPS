//
// auto-generated by ops.py on 2013-11-26 14:50
//

__constant__ int xdim0_advec_mom_kernel_x1;
__constant__ int xdim1_advec_mom_kernel_x1;
__constant__ int xdim2_advec_mom_kernel_x1;
__constant__ int xdim3_advec_mom_kernel_x1;
__constant__ int xdim4_advec_mom_kernel_x1;

#define OPS_ACC0(x,y) (x+xdim0_advec_mom_kernel_x1*(y))
#define OPS_ACC1(x,y) (x+xdim1_advec_mom_kernel_x1*(y))
#define OPS_ACC2(x,y) (x+xdim2_advec_mom_kernel_x1*(y))
#define OPS_ACC3(x,y) (x+xdim3_advec_mom_kernel_x1*(y))
#define OPS_ACC4(x,y) (x+xdim4_advec_mom_kernel_x1*(y))

//user function
__device__

inline void advec_mom_kernel_x1( double *pre_vol, double *post_vol,
                          const double *volume,
                          const double *vol_flux_x, const double *vol_flux_y) {

  post_vol[OPS_ACC1(0,0)] = volume[OPS_ACC2(0,0)] + vol_flux_y[OPS_ACC4(0,1)] -  vol_flux_y[OPS_ACC3(0,0)];
  pre_vol[OPS_ACC0(0,0)] = post_vol[OPS_ACC1(0,0)] + vol_flux_x[OPS_ACC3(1,0)] - vol_flux_x[OPS_ACC3(0,0)];

}



#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4


__global__ void ops_advec_mom_kernel_x1(
double* __restrict arg0,
double* __restrict arg1,
const double* __restrict arg2,
const double* __restrict arg3,
const double* __restrict arg4,
int size0,
int size1 ){


  int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1 + idx_y * 1 * xdim0_advec_mom_kernel_x1;
  arg1 += idx_x * 1 + idx_y * 1 * xdim1_advec_mom_kernel_x1;
  arg2 += idx_x * 1 + idx_y * 1 * xdim2_advec_mom_kernel_x1;
  arg3 += idx_x * 1 + idx_y * 1 * xdim3_advec_mom_kernel_x1;
  arg4 += idx_x * 1 + idx_y * 1 * xdim4_advec_mom_kernel_x1;

  if (idx_x < size0 && idx_y < size1) {
    advec_mom_kernel_x1(arg0, arg1, arg2, arg3,
                   arg4);
  }

}

// host stub function
void ops_par_loop_advec_mom_kernel_x1(char const *name, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4) {

  ops_arg args[5] = { arg0, arg1, arg2, arg3, arg4};


  int x_size = range[1]-range[0];
  int y_size = range[3]-range[2];

  int xdim0 = args[0].dat->block_size[0];
  int xdim1 = args[1].dat->block_size[0];
  int xdim2 = args[2].dat->block_size[0];
  int xdim3 = args[3].dat->block_size[0];
  int xdim4 = args[4].dat->block_size[0];


  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(13,"advec_mom_kernel_x1");
  ops_timers_core(&c1,&t1);

  if (OPS_kernels[13].count == 0) {
    cudaMemcpyToSymbol( xdim0_advec_mom_kernel_x1, &xdim0, sizeof(int) );
    cudaMemcpyToSymbol( xdim1_advec_mom_kernel_x1, &xdim1, sizeof(int) );
    cudaMemcpyToSymbol( xdim2_advec_mom_kernel_x1, &xdim2, sizeof(int) );
    cudaMemcpyToSymbol( xdim3_advec_mom_kernel_x1, &xdim3, sizeof(int) );
    cudaMemcpyToSymbol( xdim4_advec_mom_kernel_x1, &xdim4, sizeof(int) );
  }



  dim3 grid( (x_size-1)/OPS_block_size_x+ 1, (y_size-1)/OPS_block_size_y + 1, 1);
  dim3 block(OPS_block_size_x,OPS_block_size_y,1);




  char *p_a[5];


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


  ops_halo_exchanges_cuda(args, 5);


  //call kernel wrapper function, passing in pointers to data
  ops_advec_mom_kernel_x1<<<grid, block >>> (  (double *)p_a[0], (double *)p_a[1],
           (double *)p_a[2], (double *)p_a[3],
           (double *)p_a[4],x_size, y_size);

  if (OPS_diags>1) cudaDeviceSynchronize();
  ops_set_dirtybit_cuda(args, 5);

  //Update kernel record
  ops_timers_core(&c2,&t2);
  OPS_kernels[13].count++;
  OPS_kernels[13].time += t2-t1;
  OPS_kernels[13].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[13].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[13].transfer += ops_compute_transfer(dim, range, &arg2);
  OPS_kernels[13].transfer += ops_compute_transfer(dim, range, &arg3);
  OPS_kernels[13].transfer += ops_compute_transfer(dim, range, &arg4);
}