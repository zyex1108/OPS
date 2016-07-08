//
// auto-generated by ops.py
//
#define OPS_ACC0(x,y,z) (n_x*0+n_y*xdim0_initialise_chunk_kernel_z*0+n_z*xdim0_initialise_chunk_kernel_z*ydim0_initialise_chunk_kernel_z*1+x+xdim0_initialise_chunk_kernel_z*(y)+xdim0_initialise_chunk_kernel_z*ydim0_initialise_chunk_kernel_z*(z))
#define OPS_ACC1(x,y,z) (n_x*0+n_y*xdim1_initialise_chunk_kernel_z*0+n_z*xdim1_initialise_chunk_kernel_z*ydim1_initialise_chunk_kernel_z*1+x+xdim1_initialise_chunk_kernel_z*(y)+xdim1_initialise_chunk_kernel_z*ydim1_initialise_chunk_kernel_z*(z))
#define OPS_ACC2(x,y,z) (n_x*0+n_y*xdim2_initialise_chunk_kernel_z*0+n_z*xdim2_initialise_chunk_kernel_z*ydim2_initialise_chunk_kernel_z*1+x+xdim2_initialise_chunk_kernel_z*(y)+xdim2_initialise_chunk_kernel_z*ydim2_initialise_chunk_kernel_z*(z))


//user function

// host stub function
void ops_par_loop_initialise_chunk_kernel_z_execute(ops_kernel_descriptor *desc) {
  ops_block block = desc->block;
  int dim = desc->dim;
  int *range = desc->range;
  ops_arg arg0 = desc->args[0];
  ops_arg arg1 = desc->args[1];
  ops_arg arg2 = desc->args[2];

  //Timing
  double t1,t2,c1,c2;

  ops_arg args[3] = { arg0, arg1, arg2};



  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,3,range,51)) return;
  #endif

  //compute locally allocated range for the sub-block
  int start[3];
  int end[3];

  for ( int n=0; n<3; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }

  #ifdef OPS_DEBUG
  ops_register_args(args, "initialise_chunk_kernel_z");
  #endif



  //set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset;
  double * __restrict__ vertexz = (double *)(args[0].data + base0);

  int base1 = args[1].dat->base_offset;
  const int * __restrict__ zz = (int *)(args[1].data + base1);

  int base2 = args[2].dat->base_offset;
  double * __restrict__ vertexdz = (double *)(args[2].data + base2);


  //initialize global variable with the dimension of dats
  int xdim0_initialise_chunk_kernel_z = args[0].dat->size[0];
  int ydim0_initialise_chunk_kernel_z = args[0].dat->size[1];
  int xdim1_initialise_chunk_kernel_z = args[1].dat->size[0];
  int ydim1_initialise_chunk_kernel_z = args[1].dat->size[1];
  int xdim2_initialise_chunk_kernel_z = args[2].dat->size[0];
  int ydim2_initialise_chunk_kernel_z = args[2].dat->size[1];

  #pragma omp parallel for collapse(2)
  for ( int n_z=start[2]; n_z<end[2]; n_z++ ){
    for ( int n_y=start[1]; n_y<end[1]; n_y++ ){
      #pragma omp simd
      for ( int n_x=start[0]; n_x<end[0]; n_x++ ){
        
  int z_min=field.z_min-2;

  double min_z, d_z;
  d_z = (grid.zmax - grid.zmin)/(double)grid.z_cells;
  min_z=grid.zmin+d_z*field.back;

  vertexz[OPS_ACC0(0,0,0)] = min_z + d_z * (zz[OPS_ACC1(0,0,0)] - z_min);
  vertexdz[OPS_ACC2(0,0,0)] = (double)d_z;

      }
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2


void ops_par_loop_initialise_chunk_kernel_z(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2) {
  ops_kernel_descriptor *desc = (ops_kernel_descriptor *)malloc(sizeof(ops_kernel_descriptor));
  desc->name = name;
  desc->block = block;
  desc->dim = dim;
  desc->index = 51;
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned) return;
  for ( int n=0; n<3; n++ ){
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
  for ( int i=0; i<6; i++ ){
    desc->range[i] = range[i];
  }
  #endif //OPS_MPI
  desc->nargs = 3;
  desc->args = (ops_arg*)malloc(3*sizeof(ops_arg));
  desc->args[0] = arg0;
  desc->args[1] = arg1;
  desc->args[2] = arg2;
  desc->function = ops_par_loop_initialise_chunk_kernel_z_execute;
  ops_enqueue_kernel(desc);
  }
