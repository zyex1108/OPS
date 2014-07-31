//
// auto-generated by ops.py on 2014-07-31 12:19
//


static bool isbuilt_accelerate_kernel = false;

void buildOpenCLKernels_accelerate_kernel(int xdim0, int ydim0,
 int xdim1, int ydim1,
 int xdim2, int ydim2,
 int xdim3, int ydim3,
 int xdim4, int ydim4,
 int xdim5, int ydim5,
 int xdim6, int ydim6,
 int xdim7, int ydim7,
 int xdim8, int ydim8,
 int xdim9, int ydim9,
 int xdim10, int ydim10,
 int xdim11, int ydim11,
 int xdim12, int ydim12,
 int xdim13, int ydim13) {
  
  if(!isbuilt_accelerate_kernel) {
    buildOpenCLKernels();
    //clSafeCall( clUnloadCompiler() );       
    cl_int ret;
    char* source_filename[1] = {"./OpenCL/accelerate_kernel.cl"};

    // Load the kernel source code into the array source_str
    FILE *fid;
    char *source_str[1];
    size_t source_size[1];

    for(int i=0; i<1; i++) {
      fid = fopen(source_filename[i], "r");
      if (!fid) {
        fprintf(stderr, "Can't open the kernel source file!\n");
        exit(1);
      }
      
      source_str[i] = (char*)malloc(4*0x1000000);
      source_size[i] = fread(source_str[i], 1, 4*0x1000000, fid);
      if(source_size[i] != 4*0x1000000) {
        if (ferror(fid)) {
          printf ("Error while reading kernel source file %s\n", source_filename[i]);
          exit(-1);
        }
        if (feof(fid))
          printf ("Kernel source file %s succesfuly read.\n", source_filename[i]);
          //printf("%s\n",source_str[i]);
      }
      fclose(fid);
    }
    
    printf("Compiling accelerate_kernel source -- start \n");

      // Create a program from the source
      OPS_opencl_core.program = clCreateProgramWithSource(OPS_opencl_core.context, 1, (const char **) &source_str, (const size_t *) &source_size, &ret);
      clSafeCall( ret );

      // Build the program
      char buildOpts[255*14];
      char* pPath = NULL;
      pPath = getenv ("OPS_INSTALL_PATH");
      if (pPath!=NULL)
        sprintf(buildOpts,"-cl-mad-enable -I%s/include -DOPS_WARPSIZE=%d  -Dxdim0_accelerate_kernel=%d -Dydim0_accelerate_kernel=%d -Dxdim1_accelerate_kernel=%d -Dydim1_accelerate_kernel=%d -Dxdim2_accelerate_kernel=%d -Dydim2_accelerate_kernel=%d -Dxdim3_accelerate_kernel=%d -Dydim3_accelerate_kernel=%d -Dxdim4_accelerate_kernel=%d -Dydim4_accelerate_kernel=%d -Dxdim5_accelerate_kernel=%d -Dydim5_accelerate_kernel=%d -Dxdim6_accelerate_kernel=%d -Dydim6_accelerate_kernel=%d -Dxdim7_accelerate_kernel=%d -Dydim7_accelerate_kernel=%d -Dxdim8_accelerate_kernel=%d -Dydim8_accelerate_kernel=%d -Dxdim9_accelerate_kernel=%d -Dydim9_accelerate_kernel=%d -Dxdim10_accelerate_kernel=%d -Dydim10_accelerate_kernel=%d -Dxdim11_accelerate_kernel=%d -Dydim11_accelerate_kernel=%d -Dxdim12_accelerate_kernel=%d -Dydim12_accelerate_kernel=%d -Dxdim13_accelerate_kernel=%d -Dydim13_accelerate_kernel=%d", pPath, 32,xdim0, ydim0,xdim1, ydim1,xdim2, ydim2,xdim3, ydim3,xdim4, ydim4,xdim5, ydim5,xdim6, ydim6,xdim7, ydim7,xdim8, ydim8,xdim9, ydim9,xdim10, ydim10,xdim11, ydim11,xdim12, ydim12,xdim13, ydim13);
      else {
        sprintf("Incorrect OPS_INSTALL_PATH %s\n",pPath);
        exit(EXIT_FAILURE);
      }
      
      ret = clBuildProgram(OPS_opencl_core.program, 1, &OPS_opencl_core.device_id, buildOpts, NULL, NULL);

      if(ret != CL_SUCCESS) {
        char* build_log;
        size_t log_size;
        clSafeCall( clGetProgramBuildInfo(OPS_opencl_core.program, OPS_opencl_core.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size) );
        build_log = (char*) malloc(log_size+1);
        clSafeCall( clGetProgramBuildInfo(OPS_opencl_core.program, OPS_opencl_core.device_id, CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL) );
        build_log[log_size] = '\0';
        fprintf(stderr, "=============== OpenCL Program Build Info ================\n\n%s", build_log);
        fprintf(stderr, "\n========================================================= \n");
        free(build_log);
        exit(EXIT_FAILURE);
      }
      printf("compiling accelerate_kernel -- done\n");

    // Create the OpenCL kernel
    OPS_opencl_core.kernel[1] = clCreateKernel(OPS_opencl_core.program, "ops_accelerate_kernel", &ret);
    clSafeCall( ret );
      
    isbuilt_accelerate_kernel = true;
  }
  
}  


// host stub function
void ops_par_loop_accelerate_kernel(char const *name, ops_block Block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6, ops_arg arg7, ops_arg arg8,
 ops_arg arg9, ops_arg arg10, ops_arg arg11, ops_arg arg12, ops_arg arg13) {
  ops_arg args[14] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13};

  sub_block_list sb = OPS_sub_block_list[Block->index];
  //compute localy allocated range for the sub-block
  int start_add[3];
  int end_add[3];
  for ( int n=0; n<3; n++ ){
    start_add[n] = sb->istart[n];end_add[n] = sb->iend[n]+1;
    if (start_add[n] >= range[2*n]) {
      start_add[n] = 0;
    }
    else {
      start_add[n] = range[2*n] - start_add[n];
    }
    if (end_add[n] >= range[2*n+1]) {
      end_add[n] = range[2*n+1] - sb->istart[n];
    }
    else {
      end_add[n] = sb->sizes[n];
    }
  }


  int x_size = MAX(0,end_add[0]-start_add[0]);
  int y_size = MAX(0,end_add[1]-start_add[1]);
  int z_size = MAX(0,end_add[2]-start_add[2]);

  int xdim0 = args[0].dat->block_size[0]*args[0].dat->dim;
  int ydim0 = args[0].dat->block_size[1];
  int xdim1 = args[1].dat->block_size[0]*args[1].dat->dim;
  int ydim1 = args[1].dat->block_size[1];
  int xdim2 = args[2].dat->block_size[0]*args[2].dat->dim;
  int ydim2 = args[2].dat->block_size[1];
  int xdim3 = args[3].dat->block_size[0]*args[3].dat->dim;
  int ydim3 = args[3].dat->block_size[1];
  int xdim4 = args[4].dat->block_size[0]*args[4].dat->dim;
  int ydim4 = args[4].dat->block_size[1];
  int xdim5 = args[5].dat->block_size[0]*args[5].dat->dim;
  int ydim5 = args[5].dat->block_size[1];
  int xdim6 = args[6].dat->block_size[0]*args[6].dat->dim;
  int ydim6 = args[6].dat->block_size[1];
  int xdim7 = args[7].dat->block_size[0]*args[7].dat->dim;
  int ydim7 = args[7].dat->block_size[1];
  int xdim8 = args[8].dat->block_size[0]*args[8].dat->dim;
  int ydim8 = args[8].dat->block_size[1];
  int xdim9 = args[9].dat->block_size[0]*args[9].dat->dim;
  int ydim9 = args[9].dat->block_size[1];
  int xdim10 = args[10].dat->block_size[0]*args[10].dat->dim;
  int ydim10 = args[10].dat->block_size[1];
  int xdim11 = args[11].dat->block_size[0]*args[11].dat->dim;
  int ydim11 = args[11].dat->block_size[1];
  int xdim12 = args[12].dat->block_size[0]*args[12].dat->dim;
  int ydim12 = args[12].dat->block_size[1];
  int xdim13 = args[13].dat->block_size[0]*args[13].dat->dim;
  int ydim13 = args[13].dat->block_size[1];

  //build opencl kernel if not already built

  buildOpenCLKernels_accelerate_kernel(
  xdim0, ydim0,
  xdim1, ydim1,
  xdim2, ydim2,
  xdim3, ydim3,
  xdim4, ydim4,
  xdim5, ydim5,
  xdim6, ydim6,
  xdim7, ydim7,
  xdim8, ydim8,
  xdim9, ydim9,
  xdim10, ydim10,
  xdim11, ydim11,
  xdim12, ydim12,
  xdim13, ydim13);

  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(1,"accelerate_kernel");
  ops_timers_core(&c2,&t2);

  //set up OpenCL thread blocks
  size_t globalWorkSize[3] = {((x_size-1)/OPS_block_size_x+ 1)*OPS_block_size_x, ((y_size-1)/OPS_block_size_y + 1)*OPS_block_size_y, z_size};
  size_t localWorkSize[3] =  {OPS_block_size_x,OPS_block_size_y,1};




  int dat0 = args[0].dat->size;
  int dat1 = args[1].dat->size;
  int dat2 = args[2].dat->size;
  int dat3 = args[3].dat->size;
  int dat4 = args[4].dat->size;
  int dat5 = args[5].dat->size;
  int dat6 = args[6].dat->size;
  int dat7 = args[7].dat->size;
  int dat8 = args[8].dat->size;
  int dat9 = args[9].dat->size;
  int dat10 = args[10].dat->size;
  int dat11 = args[11].dat->size;
  int dat12 = args[12].dat->size;
  int dat13 = args[13].dat->size;

  //set up initial pointers
  int base0 = 1 * 
  (start_add[0] * args[0].stencil->stride[0] - args[0].dat->offset[0]);
  base0 = base0 + args[0].dat->block_size[0] *
  (start_add[1] * args[0].stencil->stride[1] - args[0].dat->offset[1]);
  base0 = base0 + args[0].dat->block_size[0] *  args[0].dat->block_size[1] *
  (start_add[2] * args[0].stencil->stride[2] - args[0].dat->offset[2]);

  //set up initial pointers
  int base1 = 1 * 
  (start_add[0] * args[1].stencil->stride[0] - args[1].dat->offset[0]);
  base1 = base1 + args[1].dat->block_size[0] *
  (start_add[1] * args[1].stencil->stride[1] - args[1].dat->offset[1]);
  base1 = base1 + args[1].dat->block_size[0] *  args[1].dat->block_size[1] *
  (start_add[2] * args[1].stencil->stride[2] - args[1].dat->offset[2]);

  //set up initial pointers
  int base2 = 1 * 
  (start_add[0] * args[2].stencil->stride[0] - args[2].dat->offset[0]);
  base2 = base2 + args[2].dat->block_size[0] *
  (start_add[1] * args[2].stencil->stride[1] - args[2].dat->offset[1]);
  base2 = base2 + args[2].dat->block_size[0] *  args[2].dat->block_size[1] *
  (start_add[2] * args[2].stencil->stride[2] - args[2].dat->offset[2]);

  //set up initial pointers
  int base3 = 1 * 
  (start_add[0] * args[3].stencil->stride[0] - args[3].dat->offset[0]);
  base3 = base3 + args[3].dat->block_size[0] *
  (start_add[1] * args[3].stencil->stride[1] - args[3].dat->offset[1]);
  base3 = base3 + args[3].dat->block_size[0] *  args[3].dat->block_size[1] *
  (start_add[2] * args[3].stencil->stride[2] - args[3].dat->offset[2]);

  //set up initial pointers
  int base4 = 1 * 
  (start_add[0] * args[4].stencil->stride[0] - args[4].dat->offset[0]);
  base4 = base4 + args[4].dat->block_size[0] *
  (start_add[1] * args[4].stencil->stride[1] - args[4].dat->offset[1]);
  base4 = base4 + args[4].dat->block_size[0] *  args[4].dat->block_size[1] *
  (start_add[2] * args[4].stencil->stride[2] - args[4].dat->offset[2]);

  //set up initial pointers
  int base5 = 1 * 
  (start_add[0] * args[5].stencil->stride[0] - args[5].dat->offset[0]);
  base5 = base5 + args[5].dat->block_size[0] *
  (start_add[1] * args[5].stencil->stride[1] - args[5].dat->offset[1]);
  base5 = base5 + args[5].dat->block_size[0] *  args[5].dat->block_size[1] *
  (start_add[2] * args[5].stencil->stride[2] - args[5].dat->offset[2]);

  //set up initial pointers
  int base6 = 1 * 
  (start_add[0] * args[6].stencil->stride[0] - args[6].dat->offset[0]);
  base6 = base6 + args[6].dat->block_size[0] *
  (start_add[1] * args[6].stencil->stride[1] - args[6].dat->offset[1]);
  base6 = base6 + args[6].dat->block_size[0] *  args[6].dat->block_size[1] *
  (start_add[2] * args[6].stencil->stride[2] - args[6].dat->offset[2]);

  //set up initial pointers
  int base7 = 1 * 
  (start_add[0] * args[7].stencil->stride[0] - args[7].dat->offset[0]);
  base7 = base7 + args[7].dat->block_size[0] *
  (start_add[1] * args[7].stencil->stride[1] - args[7].dat->offset[1]);
  base7 = base7 + args[7].dat->block_size[0] *  args[7].dat->block_size[1] *
  (start_add[2] * args[7].stencil->stride[2] - args[7].dat->offset[2]);

  //set up initial pointers
  int base8 = 1 * 
  (start_add[0] * args[8].stencil->stride[0] - args[8].dat->offset[0]);
  base8 = base8 + args[8].dat->block_size[0] *
  (start_add[1] * args[8].stencil->stride[1] - args[8].dat->offset[1]);
  base8 = base8 + args[8].dat->block_size[0] *  args[8].dat->block_size[1] *
  (start_add[2] * args[8].stencil->stride[2] - args[8].dat->offset[2]);

  //set up initial pointers
  int base9 = 1 * 
  (start_add[0] * args[9].stencil->stride[0] - args[9].dat->offset[0]);
  base9 = base9 + args[9].dat->block_size[0] *
  (start_add[1] * args[9].stencil->stride[1] - args[9].dat->offset[1]);
  base9 = base9 + args[9].dat->block_size[0] *  args[9].dat->block_size[1] *
  (start_add[2] * args[9].stencil->stride[2] - args[9].dat->offset[2]);

  //set up initial pointers
  int base10 = 1 * 
  (start_add[0] * args[10].stencil->stride[0] - args[10].dat->offset[0]);
  base10 = base10 + args[10].dat->block_size[0] *
  (start_add[1] * args[10].stencil->stride[1] - args[10].dat->offset[1]);
  base10 = base10 + args[10].dat->block_size[0] *  args[10].dat->block_size[1] *
  (start_add[2] * args[10].stencil->stride[2] - args[10].dat->offset[2]);

  //set up initial pointers
  int base11 = 1 * 
  (start_add[0] * args[11].stencil->stride[0] - args[11].dat->offset[0]);
  base11 = base11 + args[11].dat->block_size[0] *
  (start_add[1] * args[11].stencil->stride[1] - args[11].dat->offset[1]);
  base11 = base11 + args[11].dat->block_size[0] *  args[11].dat->block_size[1] *
  (start_add[2] * args[11].stencil->stride[2] - args[11].dat->offset[2]);

  //set up initial pointers
  int base12 = 1 * 
  (start_add[0] * args[12].stencil->stride[0] - args[12].dat->offset[0]);
  base12 = base12 + args[12].dat->block_size[0] *
  (start_add[1] * args[12].stencil->stride[1] - args[12].dat->offset[1]);
  base12 = base12 + args[12].dat->block_size[0] *  args[12].dat->block_size[1] *
  (start_add[2] * args[12].stencil->stride[2] - args[12].dat->offset[2]);

  //set up initial pointers
  int base13 = 1 * 
  (start_add[0] * args[13].stencil->stride[0] - args[13].dat->offset[0]);
  base13 = base13 + args[13].dat->block_size[0] *
  (start_add[1] * args[13].stencil->stride[1] - args[13].dat->offset[1]);
  base13 = base13 + args[13].dat->block_size[0] *  args[13].dat->block_size[1] *
  (start_add[2] * args[13].stencil->stride[2] - args[13].dat->offset[2]);


  ops_H_D_exchanges_cuda(args, 14);

  ops_timers_core(&c1,&t1);
  OPS_kernels[1].mpi_time += t1-t2;

  clSafeCall( clEnqueueWriteBuffer(OPS_opencl_core.command_queue, OPS_opencl_core.constant[14], CL_TRUE, 0, sizeof(double)*1, (void*) &dt, 0, NULL, NULL) );
  clSafeCall( clFlush(OPS_opencl_core.command_queue) );

  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 0, sizeof(cl_mem), (void*) &arg0.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 1, sizeof(cl_mem), (void*) &arg1.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 2, sizeof(cl_mem), (void*) &arg2.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 3, sizeof(cl_mem), (void*) &arg3.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 4, sizeof(cl_mem), (void*) &arg4.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 5, sizeof(cl_mem), (void*) &arg5.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 6, sizeof(cl_mem), (void*) &arg6.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 7, sizeof(cl_mem), (void*) &arg7.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 8, sizeof(cl_mem), (void*) &arg8.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 9, sizeof(cl_mem), (void*) &arg9.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 10, sizeof(cl_mem), (void*) &arg10.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 11, sizeof(cl_mem), (void*) &arg11.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 12, sizeof(cl_mem), (void*) &arg12.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 13, sizeof(cl_mem), (void*) &arg13.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 14, sizeof(cl_mem), (void*) &OPS_opencl_core.constant[14]) );
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 15, sizeof(cl_int), (void*) &base0 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 16, sizeof(cl_int), (void*) &base1 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 17, sizeof(cl_int), (void*) &base2 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 18, sizeof(cl_int), (void*) &base3 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 19, sizeof(cl_int), (void*) &base4 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 20, sizeof(cl_int), (void*) &base5 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 21, sizeof(cl_int), (void*) &base6 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 22, sizeof(cl_int), (void*) &base7 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 23, sizeof(cl_int), (void*) &base8 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 24, sizeof(cl_int), (void*) &base9 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 25, sizeof(cl_int), (void*) &base10 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 26, sizeof(cl_int), (void*) &base11 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 27, sizeof(cl_int), (void*) &base12 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 28, sizeof(cl_int), (void*) &base13 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 29, sizeof(cl_int), (void*) &x_size ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 30, sizeof(cl_int), (void*) &y_size ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[1], 31, sizeof(cl_int), (void*) &z_size ));

  //call/enque opencl kernel wrapper function
  clSafeCall( clEnqueueNDRangeKernel(OPS_opencl_core.command_queue, OPS_opencl_core.kernel[1], 3, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL) );
  if (OPS_diags>1) {
    clSafeCall( clFinish(OPS_opencl_core.command_queue) );
  }

  ops_set_dirtybit_opencl(args, 14);

  //Update kernel record
  ops_timers_core(&c2,&t2);
  OPS_kernels[1].count++;
  OPS_kernels[1].time += t2-t1;
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg2);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg3);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg4);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg5);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg6);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg7);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg8);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg9);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg10);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg11);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg12);
  OPS_kernels[1].transfer += ops_compute_transfer(dim, range, &arg13);
}
