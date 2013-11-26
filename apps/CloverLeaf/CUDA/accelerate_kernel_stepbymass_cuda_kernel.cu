//
// auto-generated by ops.py on 2013-11-26 11:43
//

__constant__ int xdim0_accelerate_kernel_stepbymass;
__constant__ int xdim1_accelerate_kernel_stepbymass;
__constant__ int xdim2_accelerate_kernel_stepbymass;

#define OPS_ACC0(x,y) (x+xdim0_accelerate_kernel_stepbymass*(y))
#define OPS_ACC1(x,y) (x+xdim1_accelerate_kernel_stepbymass*(y))
#define OPS_ACC2(x,y) (x+xdim2_accelerate_kernel_stepbymass*(y))

//user function
__device__

inline void accelerate_kernel_stepbymass(const double *density0, const double *volume,
                double *stepbymass) {

  double nodal_mass;

  nodal_mass = ( density0[OPS_ACC0(-1,-1)] * volume[OPS_ACC1(-1,-1)]
    + density0[OPS_ACC0(0,-1)] * volume[OPS_ACC1(0,-1)]
    + density0[OPS_ACC0(0,0)] * volume[OPS_ACC1(0,0)]
    + density0[OPS_ACC0(-1,0)] * volume[OPS_ACC1(-1,0)] ) * 0.25;

  stepbymass[OPS_ACC2(0,0)] = 0.5*dt / nodal_mass;

}



#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2


__global__ void ops_accelerate_kernel_stepbymass(
const double* __restrict arg0,
const double* __restrict arg1,
double* __restrict arg2,
int size0,
int size1 ){


  int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1 + idx_y * 1 * xdim0_accelerate_kernel_stepbymass;
  arg1 += idx_x * 1 + idx_y * 1 * xdim1_accelerate_kernel_stepbymass;
  arg2 += idx_x * 1 + idx_y * 1 * xdim2_accelerate_kernel_stepbymass;

  if (idx_x < size0 && idx_y < size1) {
    accelerate_kernel_stepbymass(arg0, arg1, arg2);
  }

}

// host stub function
void ops_par_loop_accelerate_kernel_stepbymass(char const *name, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2) {

  ops_arg args[3] = { arg0, arg1, arg2};


  int x_size = range[1]-range[0];
  int y_size = range[3]-range[2];

  int xdim0 = args[0].dat->block_size[0];
  int xdim1 = args[1].dat->block_size[0];
  int xdim2 = args[2].dat->block_size[0];


  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(28,"accelerate_kernel_stepbymass");
  ops_timers_core(&c1,&t1);

  if (OPS_kernels[28].count == 0) {
    cudaMemcpyToSymbol( xdim0_accelerate_kernel_stepbymass, &xdim0, sizeof(int) );
    cudaMemcpyToSymbol( xdim1_accelerate_kernel_stepbymass, &xdim1, sizeof(int) );
    cudaMemcpyToSymbol( xdim2_accelerate_kernel_stepbymass, &xdim2, sizeof(int) );
  }



  dim3 grid( (x_size-1)/OPS_block_size_x+ 1, (y_size-1)/OPS_block_size_y + 1, 1);
  dim3 block(OPS_block_size_x,OPS_block_size_y,1);




  char *p_a[3];


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


  ops_halo_exchanges_cuda(args, 3);


  //call kernel wrapper function, passing in pointers to data
  ops_accelerate_kernel_stepbymass<<<grid, block >>> (  (double *)p_a[0], (double *)p_a[1],
           (double *)p_a[2],x_size, y_size);

  if (OPS_diags>1) cudaDeviceSynchronize();
  ops_set_dirtybit_cuda(args, 3);

  //Update kernel record
  ops_timers_core(&c2,&t2);
  OPS_kernels[28].count++;
  OPS_kernels[28].time += t2-t1;
  OPS_kernels[28].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[28].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[28].transfer += ops_compute_transfer(dim, range, &arg2);
}
