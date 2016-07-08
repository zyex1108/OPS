//
// auto-generated by ops.py
//
#include "./OpenACC/clover_leaf_common.h"

#undef OPS_GPU

extern int xdim0_initialise_chunk_kernel_volume;
int xdim0_initialise_chunk_kernel_volume_h = -1;
extern int ydim0_initialise_chunk_kernel_volume;
int ydim0_initialise_chunk_kernel_volume_h = -1;
extern int xdim1_initialise_chunk_kernel_volume;
int xdim1_initialise_chunk_kernel_volume_h = -1;
extern int ydim1_initialise_chunk_kernel_volume;
int ydim1_initialise_chunk_kernel_volume_h = -1;
extern int xdim2_initialise_chunk_kernel_volume;
int xdim2_initialise_chunk_kernel_volume_h = -1;
extern int ydim2_initialise_chunk_kernel_volume;
int ydim2_initialise_chunk_kernel_volume_h = -1;
extern int xdim3_initialise_chunk_kernel_volume;
int xdim3_initialise_chunk_kernel_volume_h = -1;
extern int ydim3_initialise_chunk_kernel_volume;
int ydim3_initialise_chunk_kernel_volume_h = -1;
extern int xdim4_initialise_chunk_kernel_volume;
int xdim4_initialise_chunk_kernel_volume_h = -1;
extern int ydim4_initialise_chunk_kernel_volume;
int ydim4_initialise_chunk_kernel_volume_h = -1;
extern int xdim5_initialise_chunk_kernel_volume;
int xdim5_initialise_chunk_kernel_volume_h = -1;
extern int ydim5_initialise_chunk_kernel_volume;
int ydim5_initialise_chunk_kernel_volume_h = -1;
extern int xdim6_initialise_chunk_kernel_volume;
int xdim6_initialise_chunk_kernel_volume_h = -1;
extern int ydim6_initialise_chunk_kernel_volume;
int ydim6_initialise_chunk_kernel_volume_h = -1;

#ifdef __cplusplus
extern "C" {
#endif
void initialise_chunk_kernel_volume_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  double *p_a3,
  double *p_a4,
  double *p_a5,
  double *p_a6,
  int x_size, int y_size, int z_size);

#ifdef __cplusplus
}
#endif

// host stub function
void ops_par_loop_initialise_chunk_kernel_volume(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3, ops_arg arg4, ops_arg arg5, ops_arg arg6) {

  //Timing
  double t1,t2,c1,c2;
  ops_arg args[7] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6};


  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,7,range,55)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(55,"initialise_chunk_kernel_volume");
    OPS_kernels[55].count++;
    ops_timers_core(&c1,&t1);
  }

  //compute localy allocated range for the sub-block

  int start[3];
  int end[3];
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned) return;
  for ( int n=0; n<3; n++ ){
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
  for ( int n=0; n<3; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }
  #endif //OPS_MPI

  int x_size = MAX(0,end[0]-start[0]);
  int y_size = MAX(0,end[1]-start[1]);
  int z_size = MAX(0,end[2]-start[2]);


  xdim0 = args[0].dat->size[0];
  ydim0 = args[0].dat->size[1];
  xdim1 = args[1].dat->size[0];
  ydim1 = args[1].dat->size[1];
  xdim2 = args[2].dat->size[0];
  ydim2 = args[2].dat->size[1];
  xdim3 = args[3].dat->size[0];
  ydim3 = args[3].dat->size[1];
  xdim4 = args[4].dat->size[0];
  ydim4 = args[4].dat->size[1];
  xdim5 = args[5].dat->size[0];
  ydim5 = args[5].dat->size[1];
  xdim6 = args[6].dat->size[0];
  ydim6 = args[6].dat->size[1];
  if (xdim0 != xdim0_initialise_chunk_kernel_volume_h || ydim0 != ydim0_initialise_chunk_kernel_volume_h || xdim1 != xdim1_initialise_chunk_kernel_volume_h || ydim1 != ydim1_initialise_chunk_kernel_volume_h || xdim2 != xdim2_initialise_chunk_kernel_volume_h || ydim2 != ydim2_initialise_chunk_kernel_volume_h || xdim3 != xdim3_initialise_chunk_kernel_volume_h || ydim3 != ydim3_initialise_chunk_kernel_volume_h || xdim4 != xdim4_initialise_chunk_kernel_volume_h || ydim4 != ydim4_initialise_chunk_kernel_volume_h || xdim5 != xdim5_initialise_chunk_kernel_volume_h || ydim5 != ydim5_initialise_chunk_kernel_volume_h || xdim6 != xdim6_initialise_chunk_kernel_volume_h || ydim6 != ydim6_initialise_chunk_kernel_volume_h) {
    xdim0_initialise_chunk_kernel_volume = xdim0;
    xdim0_initialise_chunk_kernel_volume_h = xdim0;
    ydim0_initialise_chunk_kernel_volume = ydim0;
    ydim0_initialise_chunk_kernel_volume_h = ydim0;
    xdim1_initialise_chunk_kernel_volume = xdim1;
    xdim1_initialise_chunk_kernel_volume_h = xdim1;
    ydim1_initialise_chunk_kernel_volume = ydim1;
    ydim1_initialise_chunk_kernel_volume_h = ydim1;
    xdim2_initialise_chunk_kernel_volume = xdim2;
    xdim2_initialise_chunk_kernel_volume_h = xdim2;
    ydim2_initialise_chunk_kernel_volume = ydim2;
    ydim2_initialise_chunk_kernel_volume_h = ydim2;
    xdim3_initialise_chunk_kernel_volume = xdim3;
    xdim3_initialise_chunk_kernel_volume_h = xdim3;
    ydim3_initialise_chunk_kernel_volume = ydim3;
    ydim3_initialise_chunk_kernel_volume_h = ydim3;
    xdim4_initialise_chunk_kernel_volume = xdim4;
    xdim4_initialise_chunk_kernel_volume_h = xdim4;
    ydim4_initialise_chunk_kernel_volume = ydim4;
    ydim4_initialise_chunk_kernel_volume_h = ydim4;
    xdim5_initialise_chunk_kernel_volume = xdim5;
    xdim5_initialise_chunk_kernel_volume_h = xdim5;
    ydim5_initialise_chunk_kernel_volume = ydim5;
    ydim5_initialise_chunk_kernel_volume_h = ydim5;
    xdim6_initialise_chunk_kernel_volume = xdim6;
    xdim6_initialise_chunk_kernel_volume_h = xdim6;
    ydim6_initialise_chunk_kernel_volume = ydim6;
    ydim6_initialise_chunk_kernel_volume_h = ydim6;
  }



  //set up initial pointers
  int base0 = args[0].dat->base_offset + args[0].dat->elem_size * start[0] * args[0].stencil->stride[0];
  base0 = base0 + args[0].dat->elem_size *
    args[0].dat->size[0] *
    start[1] * args[0].stencil->stride[1];
  base0 = base0 + args[0].dat->elem_size *
    args[0].dat->size[0] *
    args[0].dat->size[1] *
    start[2] * args[0].stencil->stride[2];
  #ifdef OPS_GPU
  double *p_a0 = (double *)((char *)args[0].data_d + base0);
  #else
  double *p_a0 = (double *)((char *)args[0].data + base0);
  #endif

  int base1 = args[1].dat->base_offset + args[1].dat->elem_size * start[0] * args[1].stencil->stride[0];
  base1 = base1 + args[1].dat->elem_size *
    args[1].dat->size[0] *
    start[1] * args[1].stencil->stride[1];
  base1 = base1 + args[1].dat->elem_size *
    args[1].dat->size[0] *
    args[1].dat->size[1] *
    start[2] * args[1].stencil->stride[2];
  #ifdef OPS_GPU
  double *p_a1 = (double *)((char *)args[1].data_d + base1);
  #else
  double *p_a1 = (double *)((char *)args[1].data + base1);
  #endif

  int base2 = args[2].dat->base_offset + args[2].dat->elem_size * start[0] * args[2].stencil->stride[0];
  base2 = base2 + args[2].dat->elem_size *
    args[2].dat->size[0] *
    start[1] * args[2].stencil->stride[1];
  base2 = base2 + args[2].dat->elem_size *
    args[2].dat->size[0] *
    args[2].dat->size[1] *
    start[2] * args[2].stencil->stride[2];
  #ifdef OPS_GPU
  double *p_a2 = (double *)((char *)args[2].data_d + base2);
  #else
  double *p_a2 = (double *)((char *)args[2].data + base2);
  #endif

  int base3 = args[3].dat->base_offset + args[3].dat->elem_size * start[0] * args[3].stencil->stride[0];
  base3 = base3 + args[3].dat->elem_size *
    args[3].dat->size[0] *
    start[1] * args[3].stencil->stride[1];
  base3 = base3 + args[3].dat->elem_size *
    args[3].dat->size[0] *
    args[3].dat->size[1] *
    start[2] * args[3].stencil->stride[2];
  #ifdef OPS_GPU
  double *p_a3 = (double *)((char *)args[3].data_d + base3);
  #else
  double *p_a3 = (double *)((char *)args[3].data + base3);
  #endif

  int base4 = args[4].dat->base_offset + args[4].dat->elem_size * start[0] * args[4].stencil->stride[0];
  base4 = base4 + args[4].dat->elem_size *
    args[4].dat->size[0] *
    start[1] * args[4].stencil->stride[1];
  base4 = base4 + args[4].dat->elem_size *
    args[4].dat->size[0] *
    args[4].dat->size[1] *
    start[2] * args[4].stencil->stride[2];
  #ifdef OPS_GPU
  double *p_a4 = (double *)((char *)args[4].data_d + base4);
  #else
  double *p_a4 = (double *)((char *)args[4].data + base4);
  #endif

  int base5 = args[5].dat->base_offset + args[5].dat->elem_size * start[0] * args[5].stencil->stride[0];
  base5 = base5 + args[5].dat->elem_size *
    args[5].dat->size[0] *
    start[1] * args[5].stencil->stride[1];
  base5 = base5 + args[5].dat->elem_size *
    args[5].dat->size[0] *
    args[5].dat->size[1] *
    start[2] * args[5].stencil->stride[2];
  #ifdef OPS_GPU
  double *p_a5 = (double *)((char *)args[5].data_d + base5);
  #else
  double *p_a5 = (double *)((char *)args[5].data + base5);
  #endif

  int base6 = args[6].dat->base_offset + args[6].dat->elem_size * start[0] * args[6].stencil->stride[0];
  base6 = base6 + args[6].dat->elem_size *
    args[6].dat->size[0] *
    start[1] * args[6].stencil->stride[1];
  base6 = base6 + args[6].dat->elem_size *
    args[6].dat->size[0] *
    args[6].dat->size[1] *
    start[2] * args[6].stencil->stride[2];
  #ifdef OPS_GPU
  double *p_a6 = (double *)((char *)args[6].data_d + base6);
  #else
  double *p_a6 = (double *)((char *)args[6].data + base6);
  #endif


  #ifdef OPS_GPU
  ops_H_D_exchanges_device(args, 7);
  #else
  ops_H_D_exchanges_host(args, 7);
  #endif
  ops_halo_exchanges(args,7,range);

  #ifdef OPS_GPU
  ops_H_D_exchanges_device(args, 7);
  #else
  ops_H_D_exchanges_host(args, 7);
  #endif
  if (OPS_diags > 1) {
    ops_timers_core(&c2,&t2);
    OPS_kernels[55].mpi_time += t2-t1;
  }

  initialise_chunk_kernel_volume_c_wrapper(
    p_a0,
    p_a1,
    p_a2,
    p_a3,
    p_a4,
    p_a5,
    p_a6,
    x_size, y_size, z_size);

  if (OPS_diags > 1) {
    ops_timers_core(&c1,&t1);
    OPS_kernels[55].time += t1-t2;
  }
  #ifdef OPS_GPU
  ops_set_dirtybit_device(args, 7);
  #else
  ops_set_dirtybit_host(args, 7);
  #endif
  ops_set_halo_dirtybit3(&args[0],range);
  ops_set_halo_dirtybit3(&args[2],range);
  ops_set_halo_dirtybit3(&args[4],range);
  ops_set_halo_dirtybit3(&args[6],range);

  if (OPS_diags > 1) {
    //Update kernel record
    ops_timers_core(&c2,&t2);
    OPS_kernels[55].mpi_time += t2-t1;
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg0);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg1);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg2);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg3);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg4);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg5);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg6);
  }
}
