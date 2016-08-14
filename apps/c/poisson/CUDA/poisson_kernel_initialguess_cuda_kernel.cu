//
// auto-generated by ops.py
//
__constant__ int xdim0_poisson_kernel_initialguess;
int xdim0_poisson_kernel_initialguess_h = -1;
int ydim0_poisson_kernel_initialguess_h = -1;

#undef OPS_ACC0


#define OPS_ACC0(x,y) (x+xdim0_poisson_kernel_initialguess*(y))

//user function
__device__

void poisson_kernel_initialguess(double *u) {
  u[OPS_ACC0(0,0)] = 0.0;
}



#undef OPS_ACC0


__global__ void ops_poisson_kernel_initialguess(
double* __restrict arg0,
int size0,
int size1 ){


  int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1*1 + idx_y * 1*1 * xdim0_poisson_kernel_initialguess;

  if (idx_x < size0 && idx_y < size1) {
    poisson_kernel_initialguess(arg0);
  }

}

// host stub function
void ops_par_loop_poisson_kernel_initialguess(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0) {

  //Timing
  double t1,t2,c1,c2;

  ops_arg args[1] = { arg0};


  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,1,range,2)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(2,"poisson_kernel_initialguess");
    OPS_kernels[2].count++;
    ops_timers_core(&c1,&t1);
  }

  //compute locally allocated range for the sub-block
  int start[2];
  int end[2];
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned) return;
  for ( int n=0; n<2; n++ ){
    start[n] = sb->decomp_disp[n];end[n] = sb->decomp_disp[n]+sb->decomp_size[n];
    if (start[n] >= range[2*n]) {
      start[n] = 0;
    }
    else {
      start[n] = range[2*n] - start[n];
    }
    if (sb->id_m[n]==MPI_PROC_NULL && range[2*n] < 0) start[n] = range[2*n];
    if (end[n] >= range[2*n+1]) {
      end[n] = range[2*n+1] - sb->decomp_disp[n];
    }
    else {
      end[n] = sb->decomp_size[n];
    }
    if (sb->id_p[n]==MPI_PROC_NULL && (range[2*n+1] > sb->decomp_disp[n]+sb->decomp_size[n]))
      end[n] += (range[2*n+1]-sb->decomp_disp[n]-sb->decomp_size[n]);
  }
  #else //OPS_MPI
  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }
  #endif //OPS_MPI

  int x_size = MAX(0,end[0]-start[0]);
  int y_size = MAX(0,end[1]-start[1]);

  int xdim0 = args[0].dat->size[0];

  if (xdim0 != xdim0_poisson_kernel_initialguess_h) {
    cudaMemcpyToSymbol( xdim0_poisson_kernel_initialguess, &xdim0, sizeof(int) );
    xdim0_poisson_kernel_initialguess_h = xdim0;
  }



  dim3 grid( (x_size-1)/OPS_block_size_x+ 1, (y_size-1)/OPS_block_size_y + 1, 1);
  dim3 tblock(OPS_block_size_x,OPS_block_size_y,1);



  int dat0 = args[0].dat->elem_size;

  char *p_a[1];

  //set up initial pointers
  int d_m[OPS_MAX_DIM];
  #ifdef OPS_MPI
  for (int d = 0; d < dim; d++) d_m[d] = args[0].dat->d_m[d] + OPS_sub_dat_list[args[0].dat->index]->d_im[d];
  #else //OPS_MPI
  for (int d = 0; d < dim; d++) d_m[d] = args[0].dat->d_m[d];
  #endif //OPS_MPI
  int base0 = dat0 * 1 *
  (start[0] * args[0].stencil->stride[0] - args[0].dat->base[0] - d_m[0]);
  base0 = base0+ dat0 *
    args[0].dat->size[0] *
    (start[1] * args[0].stencil->stride[1] - args[0].dat->base[1] - d_m[1]);
  p_a[0] = (char *)args[0].data_d + base0;


  ops_H_D_exchanges_device(args, 1);
  ops_halo_exchanges(args,1,range);

  if (OPS_diags > 1) {
    ops_timers_core(&c2,&t2);
    OPS_kernels[2].mpi_time += t2-t1;
  }


  //call kernel wrapper function, passing in pointers to data
  ops_poisson_kernel_initialguess<<<grid, tblock >>> (  (double *)p_a[0],x_size, y_size);

  if (OPS_diags>1) {
    cutilSafeCall(cudaDeviceSynchronize());
    ops_timers_core(&c1,&t1);
    OPS_kernels[2].time += t1-t2;
  }

  ops_set_dirtybit_device(args, 1);
  ops_set_halo_dirtybit3(&args[0],range);

  if (OPS_diags > 1) {
    //Update kernel record
    ops_timers_core(&c2,&t2);
    OPS_kernels[2].mpi_time += t2-t1;
    OPS_kernels[2].transfer += ops_compute_transfer(dim, start, end, &arg0);
  }
}
