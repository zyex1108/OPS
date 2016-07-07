//
// auto-generated by ops.py
//

//user function
inline void advec_mom_kernel1_y_nonvector( const double *node_flux, const double *node_mass_pre,
                       double *mom_flux,
                       const double *celldy, const double *vel1) {





  double sigma, wind, width;
  double vdiffuw, vdiffdw, auw, adw, limiter;
  int upwind, donor, downwind, dif;
  double advec_vel_temp;

  if( (node_flux[OPS_ACC0(0,0)]) < 0.0) {
    upwind = 2;
    donor = 1;
    downwind = 0;
    dif = donor;
  } else {
    upwind = -1;
    donor = 0;
    downwind = 1;
    dif = upwind;
  }

  sigma = fabs(node_flux[OPS_ACC0(0,0)])/node_mass_pre[OPS_ACC1(0,donor)];
  width = celldy[OPS_ACC3(0,0)];
  vdiffuw = vel1[OPS_ACC4(0,donor)] - vel1[OPS_ACC4(0,upwind)];
  vdiffdw = vel1[OPS_ACC4(0,downwind)] - vel1[OPS_ACC4(0,donor)];
  limiter = 0.0;
  if(vdiffuw*vdiffdw > 0.0) {
    auw = fabs(vdiffuw);
    adw = fabs(vdiffdw);
    wind = 1.0;
    if(vdiffdw <= 0.0) wind = -1.0;
    limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy[OPS_ACC3(0,dif)])/6.0,MIN(auw,adw));
  }
  advec_vel_temp= vel1[OPS_ACC4(0,donor)] + (1.0 - sigma) * limiter;
  mom_flux[OPS_ACC2(0,0)] = advec_vel_temp * node_flux[OPS_ACC0(0,0)];
}





// host stub function
void ops_par_loop_advec_mom_kernel1_y_nonvector_execute(ops_kernel_descriptor *desc) {
  ops_block block = desc->block;
  int dim = desc->dim;
  int *range = desc->range;
  ops_arg arg0 = desc->args[0];
  ops_arg arg1 = desc->args[1];
  ops_arg arg2 = desc->args[2];
  ops_arg arg3 = desc->args[3];
  ops_arg arg4 = desc->args[4];

  //Timing
  double t1,t2,c1,c2;

  char *p_a[5];
  int  offs[5][2];
  ops_arg args[5] = { arg0, arg1, arg2, arg3, arg4};



  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,5,range,25)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(25,"advec_mom_kernel1_y_nonvector");
    OPS_kernels[25].count++;
    ops_timers_core(&c2,&t2);
  }

  //compute locally allocated range for the sub-block
  int start[2];
  int end[2];

  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }

  #ifdef OPS_DEBUG
  ops_register_args(args, "advec_mom_kernel1_y_nonvector");
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


  //initialize global variable with the dimension of dats
  xdim0 = args[0].dat->size[0];
  xdim1 = args[1].dat->size[0];
  xdim2 = args[2].dat->size[0];
  xdim3 = args[3].dat->size[0];
  xdim4 = args[4].dat->size[0];

  #pragma omp parallel for
  for ( int n_y=start[1]; n_y<end[1]; n_y++ ){
    #pragma omp simd
    for ( int n_x=start[0]; n_x<end[0]; n_x++ ){
      advec_mom_kernel1_y_nonvector(  p_a0 + n_x*1*1 + n_y*xdim0*1*1,
           p_a1 + n_x*1*1 + n_y*xdim1*1*1, p_a2 + n_x*1*1 + n_y*xdim2*1*1, p_a3 + n_x*0*1 + n_y*xdim3*1*1, p_a4 + n_x*1*1 + n_y*xdim4*1*1 );

    }
  }
}

void ops_par_loop_advec_mom_kernel1_y_nonvector(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4) {
  ops_kernel_descriptor *desc = (ops_kernel_descriptor *)malloc(sizeof(ops_kernel_descriptor));
  desc->name = name;
  desc->block = block;
  desc->dim = dim;
  desc->index = 25;
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
  desc->nargs = 5;
  desc->args = (ops_arg*)malloc(5*sizeof(ops_arg));
  desc->args[0] = arg0;
  desc->args[1] = arg1;
  desc->args[2] = arg2;
  desc->args[3] = arg3;
  desc->args[4] = arg4;
  desc->function = ops_par_loop_advec_mom_kernel1_y_nonvector_execute;
  ops_enqueue_kernel(desc);
  }
