//
// auto-generated by ops.py
//

//user function
inline void update_halo_kernel1_t2(double *density0, double *density1,
                          double *energy0, double *energy1,
                          double *pressure, double *viscosity,
                          double *soundspeed , const int* fields) {
  if(fields[FIELD_DENSITY0] == 1) density0[OPS_ACC0(0,0)] = density0[OPS_ACC0(0,-3)];
  if(fields[FIELD_DENSITY1] == 1) density1[OPS_ACC1(0,0)] = density1[OPS_ACC1(0,-3)];
  if(fields[FIELD_ENERGY0] == 1) energy0[OPS_ACC2(0,0)] = energy0[OPS_ACC2(0,-3)];
  if(fields[FIELD_ENERGY1] == 1) energy1[OPS_ACC3(0,0)] = energy1[OPS_ACC3(0,-3)];
  if(fields[FIELD_PRESSURE] == 1) pressure[OPS_ACC4(0,0)] = pressure[OPS_ACC4(0,-3)];
  if(fields[FIELD_VISCOSITY] == 1) viscosity[OPS_ACC5(0,0)] = viscosity[OPS_ACC5(0,-3)];
  if(fields[FIELD_SOUNDSPEED] == 1) soundspeed[OPS_ACC6(0,0)] = soundspeed[OPS_ACC6(0,-3)];

}





// host stub function
void ops_par_loop_update_halo_kernel1_t2_execute(ops_kernel_descriptor *desc) {
  ops_block block = desc->block;
  int dim = desc->dim;
  int *range = desc->range;
  ops_arg arg0 = desc->args[0];
  ops_arg arg1 = desc->args[1];
  ops_arg arg2 = desc->args[2];
  ops_arg arg3 = desc->args[3];
  ops_arg arg4 = desc->args[4];
  ops_arg arg5 = desc->args[5];
  ops_arg arg6 = desc->args[6];
  ops_arg arg7 = desc->args[7];

  //Timing
  double t1,t2,c1,c2;

  char *p_a[8];
  int  offs[8][2];
  ops_arg args[8] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7};



  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,8,range,45)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(45,"update_halo_kernel1_t2");
    OPS_kernels[45].count++;
    ops_timers_core(&c2,&t2);
  }

  //compute locally allocated range for the sub-block
  int start[2];
  int end[2];

  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }

  #ifdef OPS_DEBUG
  ops_register_args(args, "update_halo_kernel1_t2");
  #endif



  //set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset;
  double * __restrict__ p_a0 = (double *)(args[0].data + base0);

  int base1 = args[1].dat->base_offset;
  double * __restrict__ p_a1 = (double *)(args[1].data + base1);

  int base2 = args[2].dat->base_offset;
  double * __restrict__ p_a2 = (double *)(args[2].data + base2);

  int base3 = args[3].dat->base_offset;
  double * __restrict__ p_a3 = (double *)(args[3].data + base3);

  int base4 = args[4].dat->base_offset;
  double * __restrict__ p_a4 = (double *)(args[4].data + base4);

  int base5 = args[5].dat->base_offset;
  double * __restrict__ p_a5 = (double *)(args[5].data + base5);

  int base6 = args[6].dat->base_offset;
  double * __restrict__ p_a6 = (double *)(args[6].data + base6);

  int * __restrict__ p_a7 = (int *)args[7].data;



  //initialize global variable with the dimension of dats
  xdim0 = args[0].dat->size[0];
  xdim1 = args[1].dat->size[0];
  xdim2 = args[2].dat->size[0];
  xdim3 = args[3].dat->size[0];
  xdim4 = args[4].dat->size[0];
  xdim5 = args[5].dat->size[0];
  xdim6 = args[6].dat->size[0];

  #pragma omp parallel for
  for ( int n_y=start[1]; n_y<end[1]; n_y++ ){
    #pragma omp simd
    for ( int n_x=start[0]; n_x<end[0]; n_x++ ){
      update_halo_kernel1_t2(  p_a0 + n_x*1*1 + n_y*xdim0*1*1,
           p_a1 + n_x*1*1 + n_y*xdim1*1*1, p_a2 + n_x*1*1 + n_y*xdim2*1*1, p_a3 + n_x*1*1 + n_y*xdim3*1*1, p_a4 + n_x*1*1 + n_y*xdim4*1*1,
           p_a5 + n_x*1*1 + n_y*xdim5*1*1, p_a6 + n_x*1*1 + n_y*xdim6*1*1, p_a7 );

    }
  }
}

void ops_par_loop_update_halo_kernel1_t2(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6, ops_arg arg7) {
  ops_kernel_descriptor *desc = (ops_kernel_descriptor *)malloc(sizeof(ops_kernel_descriptor));
  desc->name = name;
  desc->block = block;
  desc->dim = dim;
  desc->index = 45;
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned) return;
  for ( int n=0; n<2; n++ ){
    desc->range[2*n] = sb->decomp_disp[n];desc->range[2*n+1] = sb->decomp_disp[n]+sb->decomp_size[n];
    if (desc->range[2*n] >= range[2*n]) {
      desc->range[2*n] = 0;
    }
    else {
      desc->range[2*n] = range[2*n] - desc->range[2*n];
    }
    if (sb->id_m[n]==MPI_PROC_NULL && range[2*n] < 0) desc->range[2*n] = range[2*n];
    if (desc->range[2*n+1] >= range[2*n+1]) {
      desc->range[2*n+1] = range[2*n+1] - sb->decomp_disp[n];
    }
    else {
      desc->range[2*n+1] = sb->decomp_size[n];
    }
    if (sb->id_p[n]==MPI_PROC_NULL && (range[2*n+1] > sb->decomp_disp[n]+sb->decomp_size[n]))
      desc->range[2*n+1] += (range[2*n+1]-sb->decomp_disp[n]-sb->decomp_size[n]);
  }
  #else //OPS_MPI
  for ( int i=0; i<4; i++ ){
    desc->range[i] = range[i];
  }
  #endif //OPS_MPI
  desc->nargs = 8;
  desc->args = (ops_arg*)malloc(8*sizeof(ops_arg));
  desc->args[0] = arg0;
  desc->args[1] = arg1;
  desc->args[2] = arg2;
  desc->args[3] = arg3;
  desc->args[4] = arg4;
  desc->args[5] = arg5;
  desc->args[6] = arg6;
  desc->args[7] = arg7;
  char *tmp = (char*)malloc(NUM_FIELDS*sizeof(int));
  memcpy(tmp, arg7.data,NUM_FIELDS*sizeof(int));
  desc->args[7].data = tmp;
  desc->function = ops_par_loop_update_halo_kernel1_t2_execute;
  ops_enqueue_kernel(desc);
  }
