//
// auto-generated by ops.py
//
__constant__ int xdim0_Riemann_kernel;
int xdim0_Riemann_kernel_h = -1;
int ydim0_Riemann_kernel_h = -1;
__constant__ int xdim1_Riemann_kernel;
int xdim1_Riemann_kernel_h = -1;
int ydim1_Riemann_kernel_h = -1;
__constant__ int xdim2_Riemann_kernel;
int xdim2_Riemann_kernel_h = -1;
int ydim2_Riemann_kernel_h = -1;
__constant__ int xdim3_Riemann_kernel;
int xdim3_Riemann_kernel_h = -1;
int ydim3_Riemann_kernel_h = -1;
__constant__ int xdim4_Riemann_kernel;
int xdim4_Riemann_kernel_h = -1;
int ydim4_Riemann_kernel_h = -1;
__constant__ int xdim5_Riemann_kernel;
int xdim5_Riemann_kernel_h = -1;
int ydim5_Riemann_kernel_h = -1;

#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2

#undef OPS_ACC_MD3
#undef OPS_ACC_MD4
#undef OPS_ACC_MD5

#define OPS_ACC0(x) (x)
#define OPS_ACC1(x) (x)
#define OPS_ACC2(x) (x)

#define OPS_ACC_MD3(d, x) ((x)*3 + (d))
#define OPS_ACC_MD4(d, x) ((x)*9 + (d))
#define OPS_ACC_MD5(d, x) ((x)*3 + (d))
// user function
__device__

    void
    Riemann_kernel(const double *rho_new, const double *rhou_new,
                   const double *rhoE_new, double *alam, double *r,
                   double *al) {

  double rl, rr, rho, u, hl, hr, h, Vsq, csq, c;
  double dw1, dw2, dw3, delpc2, rdeluc;

  rl = sqrt(rho_new[OPS_ACC0(0)]);
  rr = sqrt(rho_new[OPS_ACC0(1)]);
  rho = rl + rr;
  u = ((rhou_new[OPS_ACC1(0)] / rl) + (rhou_new[OPS_ACC1(1)] / rr)) / rho;
  double fni =
      rhou_new[OPS_ACC1(0)] * rhou_new[OPS_ACC1(0)] / rho_new[OPS_ACC0(0)];
  double p = gam1 * (rhoE_new[OPS_ACC2(0)] - 0.5 * fni);
  hl = (rhoE_new[OPS_ACC2(0)] + p) / rl;
  fni = rhou_new[OPS_ACC1(1)] * rhou_new[OPS_ACC1(1)] / rho_new[OPS_ACC0(1)];
  p = gam1 * (rhoE_new[OPS_ACC2(1)] - 0.5 * fni);
  hr = (rhoE_new[OPS_ACC2(1)] + p) / rr;
  h = (hl + hr) / rho;
  Vsq = u * u;
  csq = gam1 * (h - 0.5 * Vsq);
  c = sqrt(csq);

  alam[OPS_ACC_MD3(0, 0)] = u - c;
  alam[OPS_ACC_MD3(1, 0)] = u;
  alam[OPS_ACC_MD3(2, 0)] = u + c;

  r[OPS_ACC_MD4(0, 0)] = 1.0;
  r[OPS_ACC_MD4(1, 0)] = 1.0;
  r[OPS_ACC_MD4(2, 0)] = 1.0;

  r[OPS_ACC_MD4(3, 0)] = u - c;
  r[OPS_ACC_MD4(4, 0)] = u;
  r[OPS_ACC_MD4(5, 0)] = u + c;

  r[OPS_ACC_MD4(6, 0)] = h - u * c;
  r[OPS_ACC_MD4(7, 0)] = 0.5 * Vsq;
  r[OPS_ACC_MD4(8, 0)] = h + u * c;

  for (int m = 0; m < 9; m++)
    r[OPS_ACC_MD4(m, 0)] = r[OPS_ACC_MD4(m, 0)] / csq;

  dw1 = rho_new[OPS_ACC0(1)] - rho_new[OPS_ACC0(0)];
  dw2 = rhou_new[OPS_ACC1(1)] - rhou_new[OPS_ACC1(0)];
  dw3 = rhoE_new[OPS_ACC2(1)] - rhoE_new[OPS_ACC2(0)];

  delpc2 = gam1 * (dw3 + 0.50 * Vsq * dw1 - u * dw2) / csq;
  rdeluc = (dw2 - u * dw1) / c;

  al[OPS_ACC_MD5(0, 0)] = 0.5 * (delpc2 - rdeluc);
  al[OPS_ACC_MD5(1, 0)] = dw1 - delpc2;
  al[OPS_ACC_MD5(2, 0)] = 0.5 * (delpc2 + rdeluc);

  for (int m = 0; m < 3; m++)
    al[OPS_ACC_MD5(m, 0)] = al[OPS_ACC_MD5(m, 0)] * csq;
}

#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2

#undef OPS_ACC_MD3
#undef OPS_ACC_MD4
#undef OPS_ACC_MD5

__global__ void ops_Riemann_kernel(const double *__restrict arg0,
                                   const double *__restrict arg1,
                                   const double *__restrict arg2,
                                   double *__restrict arg3,
                                   double *__restrict arg4,
                                   double *__restrict arg5, int size0) {

  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1 * 1;
  arg1 += idx_x * 1 * 1;
  arg2 += idx_x * 1 * 1;
  arg3 += idx_x * 1 * 3;
  arg4 += idx_x * 1 * 9;
  arg5 += idx_x * 1 * 3;

  if (idx_x < size0) {
    Riemann_kernel(arg0, arg1, arg2, arg3, arg4, arg5);
  }
}

// host stub function
void ops_par_loop_Riemann_kernel(char const *name, ops_block block, int dim,
                                 int *range, ops_arg arg0, ops_arg arg1,
                                 ops_arg arg2, ops_arg arg3, ops_arg arg4,
                                 ops_arg arg5) {

  // Timing
  double t1, t2, c1, c2;

  ops_arg args[6] = {arg0, arg1, arg2, arg3, arg4, arg5};

#ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args, 6, range, 7))
    return;
#endif

  if (OPS_diags > 1) {
    ops_timing_realloc(7, "Riemann_kernel");
    OPS_kernels[7].count++;
    ops_timers_core(&c1, &t1);
  }

  // compute locally allocated range for the sub-block
  int start[1];
  int end[1];
#ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned)
    return;
  for (int n = 0; n < 1; n++) {
    start[n] = sb->decomp_disp[n];
    end[n] = sb->decomp_disp[n] + sb->decomp_size[n];
    if (start[n] >= range[2 * n]) {
      start[n] = 0;
    } else {
      start[n] = range[2 * n] - start[n];
    }
    if (sb->id_m[n] == MPI_PROC_NULL && range[2 * n] < 0)
      start[n] = range[2 * n];
    if (end[n] >= range[2 * n + 1]) {
      end[n] = range[2 * n + 1] - sb->decomp_disp[n];
    } else {
      end[n] = sb->decomp_size[n];
    }
    if (sb->id_p[n] == MPI_PROC_NULL &&
        (range[2 * n + 1] > sb->decomp_disp[n] + sb->decomp_size[n]))
      end[n] += (range[2 * n + 1] - sb->decomp_disp[n] - sb->decomp_size[n]);
  }
#else
  for (int n = 0; n < 1; n++) {
    start[n] = range[2 * n];
    end[n] = range[2 * n + 1];
  }
#endif

  int x_size = MAX(0, end[0] - start[0]);

  int xdim0 = args[0].dat->size[0];
  int xdim1 = args[1].dat->size[0];
  int xdim2 = args[2].dat->size[0];
  int xdim3 = args[3].dat->size[0];
  int xdim4 = args[4].dat->size[0];
  int xdim5 = args[5].dat->size[0];

  if (xdim0 != xdim0_Riemann_kernel_h || xdim1 != xdim1_Riemann_kernel_h ||
      xdim2 != xdim2_Riemann_kernel_h || xdim3 != xdim3_Riemann_kernel_h ||
      xdim4 != xdim4_Riemann_kernel_h || xdim5 != xdim5_Riemann_kernel_h) {
    cudaMemcpyToSymbol(xdim0_Riemann_kernel, &xdim0, sizeof(int));
    xdim0_Riemann_kernel_h = xdim0;
    cudaMemcpyToSymbol(xdim1_Riemann_kernel, &xdim1, sizeof(int));
    xdim1_Riemann_kernel_h = xdim1;
    cudaMemcpyToSymbol(xdim2_Riemann_kernel, &xdim2, sizeof(int));
    xdim2_Riemann_kernel_h = xdim2;
    cudaMemcpyToSymbol(xdim3_Riemann_kernel, &xdim3, sizeof(int));
    xdim3_Riemann_kernel_h = xdim3;
    cudaMemcpyToSymbol(xdim4_Riemann_kernel, &xdim4, sizeof(int));
    xdim4_Riemann_kernel_h = xdim4;
    cudaMemcpyToSymbol(xdim5_Riemann_kernel, &xdim5, sizeof(int));
    xdim5_Riemann_kernel_h = xdim5;
  }

  dim3 grid((x_size - 1) / OPS_block_size_x + 1, 1, 1);
  dim3 tblock(OPS_block_size_x, 1, 1);

  int dat0 = args[0].dat->elem_size;
  int dat1 = args[1].dat->elem_size;
  int dat2 = args[2].dat->elem_size;
  int dat3 = args[3].dat->elem_size;
  int dat4 = args[4].dat->elem_size;
  int dat5 = args[5].dat->elem_size;

  char *p_a[6];

  // set up initial pointers
  int d_m[OPS_MAX_DIM];
#ifdef OPS_MPI
  for (int d = 0; d < dim; d++)
    d_m[d] =
        args[0].dat->d_m[d] + OPS_sub_dat_list[args[0].dat->index]->d_im[d];
#else
  for (int d = 0; d < dim; d++)
    d_m[d] = args[0].dat->d_m[d];
#endif
  int base0 = dat0 * 1 * (start[0] * args[0].stencil->stride[0] -
                          args[0].dat->base[0] - d_m[0]);
  p_a[0] = (char *)args[0].data_d + base0;

#ifdef OPS_MPI
  for (int d = 0; d < dim; d++)
    d_m[d] =
        args[1].dat->d_m[d] + OPS_sub_dat_list[args[1].dat->index]->d_im[d];
#else
  for (int d = 0; d < dim; d++)
    d_m[d] = args[1].dat->d_m[d];
#endif
  int base1 = dat1 * 1 * (start[0] * args[1].stencil->stride[0] -
                          args[1].dat->base[0] - d_m[0]);
  p_a[1] = (char *)args[1].data_d + base1;

#ifdef OPS_MPI
  for (int d = 0; d < dim; d++)
    d_m[d] =
        args[2].dat->d_m[d] + OPS_sub_dat_list[args[2].dat->index]->d_im[d];
#else
  for (int d = 0; d < dim; d++)
    d_m[d] = args[2].dat->d_m[d];
#endif
  int base2 = dat2 * 1 * (start[0] * args[2].stencil->stride[0] -
                          args[2].dat->base[0] - d_m[0]);
  p_a[2] = (char *)args[2].data_d + base2;

#ifdef OPS_MPI
  for (int d = 0; d < dim; d++)
    d_m[d] =
        args[3].dat->d_m[d] + OPS_sub_dat_list[args[3].dat->index]->d_im[d];
#else
  for (int d = 0; d < dim; d++)
    d_m[d] = args[3].dat->d_m[d];
#endif
  int base3 = dat3 * 1 * (start[0] * args[3].stencil->stride[0] -
                          args[3].dat->base[0] - d_m[0]);
  p_a[3] = (char *)args[3].data_d + base3;

#ifdef OPS_MPI
  for (int d = 0; d < dim; d++)
    d_m[d] =
        args[4].dat->d_m[d] + OPS_sub_dat_list[args[4].dat->index]->d_im[d];
#else
  for (int d = 0; d < dim; d++)
    d_m[d] = args[4].dat->d_m[d];
#endif
  int base4 = dat4 * 1 * (start[0] * args[4].stencil->stride[0] -
                          args[4].dat->base[0] - d_m[0]);
  p_a[4] = (char *)args[4].data_d + base4;

#ifdef OPS_MPI
  for (int d = 0; d < dim; d++)
    d_m[d] =
        args[5].dat->d_m[d] + OPS_sub_dat_list[args[5].dat->index]->d_im[d];
#else
  for (int d = 0; d < dim; d++)
    d_m[d] = args[5].dat->d_m[d];
#endif
  int base5 = dat5 * 1 * (start[0] * args[5].stencil->stride[0] -
                          args[5].dat->base[0] - d_m[0]);
  p_a[5] = (char *)args[5].data_d + base5;

  ops_H_D_exchanges_device(args, 6);
  ops_halo_exchanges(args, 6, range);

  if (OPS_diags > 1) {
    ops_timers_core(&c2, &t2);
    OPS_kernels[7].mpi_time += t2 - t1;
  }

  // call kernel wrapper function, passing in pointers to data
  ops_Riemann_kernel<<<grid, tblock>>>(
      (double *)p_a[0], (double *)p_a[1], (double *)p_a[2], (double *)p_a[3],
      (double *)p_a[4], (double *)p_a[5], x_size);

  if (OPS_diags > 1) {
    cutilSafeCall(cudaDeviceSynchronize());
    ops_timers_core(&c1, &t1);
    OPS_kernels[7].time += t1 - t2;
  }

  ops_set_dirtybit_device(args, 6);
  ops_set_halo_dirtybit3(&args[3], range);
  ops_set_halo_dirtybit3(&args[4], range);
  ops_set_halo_dirtybit3(&args[5], range);

  if (OPS_diags > 1) {
    // Update kernel record
    ops_timers_core(&c2, &t2);
    OPS_kernels[7].mpi_time += t2 - t1;
    OPS_kernels[7].transfer += ops_compute_transfer(dim, start, end, &arg0);
    OPS_kernels[7].transfer += ops_compute_transfer(dim, start, end, &arg1);
    OPS_kernels[7].transfer += ops_compute_transfer(dim, start, end, &arg2);
    OPS_kernels[7].transfer += ops_compute_transfer(dim, start, end, &arg3);
    OPS_kernels[7].transfer += ops_compute_transfer(dim, start, end, &arg4);
    OPS_kernels[7].transfer += ops_compute_transfer(dim, start, end, &arg5);
  }
}
