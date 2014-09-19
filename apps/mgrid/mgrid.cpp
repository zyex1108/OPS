/*
* Open source copyright declaration based on BSD open source template:
* http://www.opensource.org/licenses/bsd-license.php
*
* This file is part of the OPS distribution.
*
* Copyright (c) 2013, Mike Giles and others. Please see the AUTHORS file in
* the main source directory for a full list of copyright holders.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
* * Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
* * The name of Mike Giles may not be used to endorse or promote products
* derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/** @Test application for multi-grid functionality
  * @author Gihan Mudalige, Istvan Reguly
  */

// standard headers
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// OPS header file
#include "ops_seq_mgrid.h"

#include "mgrid_populate_kernel.h"

/******************************************************************************
* Main program
/******************************************************************************/
int main(int argc, char **argv)
{
  /**-------------------------- Initialisation --------------------------**/

  // OPS initialisation
  ops_init(argc,argv,2);

  //declare blocks
  ops_block grid0 = ops_decl_block(2, "grid0");
  

  //declare stencils
  int s2D_00[]         = {0,0};
  ops_stencil S2D_00 = ops_decl_stencil( 2, 1, s2D_00, "00");
  
  int str2[] = {2,2};
  int stp2[] = {1/2,1/2};
  ops_stencil S2D_STR2_00 = ops_decl_strided_stencil( 2, 1, s2D_00, str2, "00");
  ops_stencil S2D_STP2_00 = ops_decl_strided_stencil( 2, 1, s2D_00, stp2, "00");
  
  //need strided restrict stencil

  //declare datasets
  int d_p[2] = {2,2}; //max halo depths for the dat in the possitive direction
  int d_m[2] = {-2,-2}; //max halo depths for the dat in the negative direction
  int size0[2] = {10, 10}; //size of the dat -- should be identical to the block on which its define on
  int size1[2] = {5, 5}; //size of the dat -- should be identical to the block on which its define on
  int base[2] = {0,0};
  double* temp = NULL;

  ops_dat data0 = ops_decl_dat(grid0, 1, size0, base, d_m, d_p, 1 /*stride*/, temp, "double", "data0");
  ops_dat data1 = ops_decl_dat(grid0, 1, size1, base, d_m, d_p, 2 /*stride*/, temp, "double", "data1"); 

  ops_partition("");

  /**-------------------------- Computations --------------------------**/


  double ct0, ct1, et0, et1;
  ops_timers_core(&ct0, &et0);

  //populate
  int iter_range[] = {0,10,0,10};
  int iter_range_small[] = {0,5,0,5};
  
  //ops_par_loop(mgrid_populate_kernel, "mgrid_populate_kernel", grid0, 2, iter_range,
  //             ops_arg_dat(data0, S2D_00, "double", OPS_WRITE),
  //             ops_arg_idx());
  ops_par_loop(mgrid_populate_kernel, "mgrid_populate_kernel", grid0, 2, iter_range_small,   //operates on a smaller range
               ops_arg_dat(data1, S2D_00, "double", OPS_WRITE),
               ops_arg_idx());

  ops_par_loop(mgrid_prolong_kernel, "mgrid_prolong_kernel", grid0, 2, iter_range, 
               ops_arg_dat(data1, S2D_00, "double", OPS_READ),
               ops_arg_dat(data0, S2D_00, "double", OPS_WRITE),
               ops_arg_idx());
  
  ops_print_dat_to_txtfile(data0, "data.txt");  
    
  
  
  /*//do restriction loop
  ops_par_loop(mgrid_restrict_kernel, "mgrid_restrict_kernel", grid0, 2, iter_range_small, 
               ops_arg_dat(data0, S2D_STR2_00_M10_P10, "double", OPS_READ),
               ops_arg_dat(data1, S2D_00, "double", OPS_WRITE),
               ops_arg_idx());
  
  
  ops_print_dat_to_txtfile(data1, "data.txt");
  
  //do prolongation loop
  ops_par_loop(mgrid_prolong_kernel, "mgrid_prolong_kernel", grid0, 2, iter_range, 
               ops_arg_dat(data1, S2D_STP2_00_P10, "double", OPS_READ),
               ops_arg_dat(data0, S2D_00, "double", OPS_WRITE),
               ops_arg_idx());
  */
  


  ops_timers_core(&ct1, &et1);
  ops_timing_output();

  ops_printf("\nTotal Wall time %lf\n",et1-et0);

  ops_exit();
}
