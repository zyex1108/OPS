
# Open source copyright declaration based on BSD open source template:
# http://www.opensource.org/licenses/bsd-license.php
#
# This file is part of the OPS distribution.
#
# Copyright (c) 2013, Mike Giles and others. Please see the AUTHORS file in
# the main source directory for a full list of copyright holders.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# The name of Mike Giles may not be used to endorse or promote products
# derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
OPS MPI_seq code generator

This routine is called by ops.py which parses the input files

It produces a file xxx_seq_kernel.cpp for each kernel,
plus a master kernel file

"""

import re
import datetime
import os

def comm(line):
  global file_text, FORTRAN, CPP
  global depth
  prefix = ' '*depth
  if len(line) == 0:
    file_text +='\n'
  else:
    file_text +=prefix+'//'+line+'\n'

def code(text):
  global file_text, g_m
  global depth
  prefix = ''
  if len(text) != 0:
    prefix = ' '*depth

  file_text += prefix+text+'\n'

def FOR(i,start,finish):
  global file_text
  global depth
  code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+'++ ){')
  depth += 2

def FOR2(i,start,finish,increment):
  global file_text
  global depth
  code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+'+='+increment+' ){')
  depth += 2

def WHILE(line):
  global file_text
  global depth
  code('while ( '+ line+ ' ){')
  depth += 2

def ENDWHILE():
  global file_text
  global depth
  depth -= 2
  code('}')

def ENDFOR():
  global file_text
  global depth
  depth -= 2
  code('}')

def IF(line):
  global file_text
  global depth
  code('if ('+ line + ') {')
  depth += 2

def ELSEIF(line):
  global file_text
  global depth
  code('else if ('+ line + ') {')
  depth += 2

def ELSE():
  global file_text
  global depth
  code('else {')
  depth += 2

def ENDIF():
  global file_text
  global depth
  depth -= 2
  code('}')

def mult(text, i, n):
  text = text + '1'
  for nn in range (0, i):
    text = text + '* args['+str(n)+'].dat->size['+str(nn)+']'

  return text

def ops_gen_mpi(master, date, consts, kernels):

  global dims, stens
  global g_m, file_text, depth

  OPS_ID   = 1;  OPS_GBL   = 2;  OPS_MAP = 3;

  OPS_READ = 1;  OPS_WRITE = 2;  OPS_RW  = 3;
  OPS_INC  = 4;  OPS_MAX   = 5;  OPS_MIN = 6;

  accsstring = ['OPS_READ','OPS_WRITE','OPS_RW','OPS_INC','OPS_MAX','OPS_MIN' ]

  NDIM = 2 #the dimension of the application is hardcoded here .. need to get this dynamically

##########################################################################
#  create new kernel file
##########################################################################

  for nk in range (0,len(kernels)):
    arg_typ  = kernels[nk]['arg_type']
    name  = kernels[nk]['name']
    nargs = kernels[nk]['nargs']
    dim   = kernels[nk]['dim']
    dims  = kernels[nk]['dims']
    stens = kernels[nk]['stens']
    var   = kernels[nk]['var']
    accs  = kernels[nk]['accs']
    typs  = kernels[nk]['typs']
    
    NDIM = int(dim)
    #parse stencil to locate strided access
    stride = [1] * nargs * NDIM
    restrict = [1] * nargs 
    prolong = [1] * nargs

    if NDIM == 2:
      for n in range (0, nargs):
        if str(stens[n]).find('STRID2D_X') > 0:
          stride[NDIM*n+1] = 0
        elif str(stens[n]).find('STRID2D_Y') > 0:
          stride[NDIM*n] = 0

    if NDIM == 3:
      for n in range (0, nargs):
        if str(stens[n]).find('STRID3D_X') > 0:
          stride[NDIM*n+1] = 0
          stride[NDIM*n+2] = 0
        elif str(stens[n]).find('STRID3D_Y') > 0:
          stride[NDIM*n] = 0
          stride[NDIM*n+2] = 0
        elif str(stens[n]).find('STRID3D_Z') > 0:
          stride[NDIM*n] = 0
          stride[NDIM*n+1] = 0
    
    ### Determine if this is a MULTI_GRID LOOP with either restrict or prolong
    MULTI_GRID = 0
    for n in range (0, nargs):
      restrict[n] = 0
      prolong[n] = 0
      if str(stens[n]).find('RESTRICT') > 0:
        restrict[n] = 1
        MULTI_GRID = 1
      if str(stens[n]).find('PROLONG') > 0 :     
        prolong[n] = 1
        MULTI_GRID = 1

    reduction = 0
    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_gbl' and accs[n] <> OPS_READ:
        reduction = 1

    arg_idx = 0
    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_idx':
        arg_idx = 1

##########################################################################
#  start with seq kernel function
##########################################################################

    g_m = 0;
    file_text = ''
    depth = 0
    n_per_line = 4

    i = name.find('kernel')
    name2 = name[0:i-1]
    #print name2

    comm('user function')
    code('#include "'+name2+'_kernel.h"')
    comm('')
    comm(' host stub function')
    code('void ops_par_loop_'+name+'(char const *name, ops_block block, int dim, int* range,')
    text = ''
    for n in range (0, nargs):

      text = text +' ops_arg arg'+str(n)
      if nargs <> 1 and n != nargs-1:
        text = text +','
      else:
        text = text +') {'
      if n%n_per_line == 3 and n <> nargs-1:
         text = text +'\n'
    code(text);
    depth = 2

    code('');
    code('char *p_a['+str(nargs)+'];')
    code('int  offs['+str(nargs)+']['+dim+'];')

    #code('ops_printf("In loop \%s\\n","'+name+'");')

    text ='ops_arg args['+str(nargs)+'] = {'
    for n in range (0, nargs):
      text = text +' arg'+str(n)
      if nargs <> 1 and n != nargs-1:
        text = text +','
      else:
        text = text +'};\n\n'
      if n%n_per_line == 5 and n <> nargs-1:
        text = text +'\n                    '
    code(text);
    code('')
    code('ops_timing_realloc('+str(nk)+',"'+name+'");')
    code('OPS_kernels['+str(nk)+'].count++;')
    code('')
    comm('compute locally allocated range for the sub-block')
    code('int start['+str(NDIM)+'];')
    code('int end['+str(NDIM)+'];')
    if arg_idx == 1:
      code('int arg_idx['+str(NDIM)+'];')
    code('')

    code('#ifdef OPS_MPI')
    code('if (compute_ranges(args, block, range, start, end, arg_idx) < 0) return;')
    code('#else //OPS_MPI')
    FOR('n','0',str(NDIM))
    code('start[n] = range[2*n];end[n] = range[2*n+1];')
    ENDFOR()
    code('#endif //OPS_MPI')
    code('')
    code('#ifdef OPS_DEBUG')
    code('ops_register_args(args, "'+name+'");')
    code('#endif')
    code('')

    if MULTI_GRID:
      for n in range (0, nargs):
        if restrict[n]  == 1 :
          code('int start_'+str(n)+'[2]; int end_'+str(n)+'[2]; int stride_'+str(n)+'[2];')
          FOR('n','0',str(NDIM))
          code('stride_'+str(n)+'[n] = args['+str(n)+'].stencil->mgrid_stride[n];')
          code('start_'+str(n)+'[n]  = start[n];')
          code('end_'+str(n)+'[n]    = end[n];')
          ENDFOR()
        elif prolong[n] == 1:
          code('int start_'+str(n)+'[2]; int end_'+str(n)+'[2]; int stride_'+str(n)+'[2];')
          FOR('n','0',str(NDIM))
          code('stride_'+str(n)+'[n] = args['+str(n)+'].stencil->mgrid_stride[n];')
          code('start_'+str(n)+'[n]  = start[n]/stride_'+str(n)+'[n];')
          code('end_'+str(n)+'[n]    = end[n]/stride_'+str(n)+'[n];')
          ENDFOR()
      
    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_dat':
        code('offs['+str(n)+'][0] = args['+str(n)+'].stencil->stride[0]*1;  //unit step in x dimension')
        for d in range (1, NDIM):          
          if restrict[n]  == 1 or prolong[n] == 1:
            code('offs['+str(n)+']['+str(d)+'] = off'+str(NDIM)+'D('+str(d)+', &start_'+str(n)+'[0],')
            if d == 1:
              code('    &end_'+str(n)+'[0],args['+str(n)+'].dat->size, args['+str(n)+'].stencil->stride) - offs['+str(n)+']['+str(d-1)+'];')
            if d == 2:
              code('    &end_'+str(n)+'[0],args['+str(n)+'].dat->size, args['+str(n)+'].stencil->stride) - offs['+str(n)+']['+str(d-1)+'] - offs['+str(n)+']['+str(d-2)+'];')
          else:
            code('offs['+str(n)+']['+str(d)+'] = off'+str(NDIM)+'D('+str(d)+', &start[0],')
            if d == 1:
              code('    &end[0],args['+str(n)+'].dat->size, args['+str(n)+'].stencil->stride) - offs['+str(n)+']['+str(d-1)+'];')
            if d == 2:
              code('    &end[0],args['+str(n)+'].dat->size, args['+str(n)+'].stencil->stride) - offs['+str(n)+']['+str(d-1)+'] - offs['+str(n)+']['+str(d-2)+'];')
        code('')

    code('')
    if arg_idx:
      code('#ifdef OPS_MPI')
      for n in range (0,NDIM):
        code('int arg_idx_'+str(n)+' = arg_idx['+str(n)+'];')
      code('#else //OPS_MPI')
      for n in range (0,NDIM):
        code('int arg_idx_'+str(n)+' = start['+str(n)+'];')
      code('#endif //OPS_MPI')
    
    if MULTI_GRID:
      code('int global_idx['+str(NDIM)+'];')
      code('#ifdef OPS_MPI')
      for n in range (0,NDIM):
        code('global_idx['+str(n)+'] = sb->decomp_disp['+str(n)+']+start['+str(n)+'];')
      code('#else //OPS_MPI')
      for n in range (0,NDIM):
        code('global_idx['+str(n)+'] = start['+str(n)+'];')
      code('#endif //OPS_MPI')
    

    code('')

    comm('Timing')
    code('double t1,t2,c1,c2;')
    code('ops_timers_core(&c2,&t2);')
    code('')

    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_dat':
        for d in range (0, NDIM):
          code('int off'+str(n)+'_'+str(d)+' = offs['+str(n)+']['+str(d)+'];')
        code('int dat'+str(n)+' = args['+str(n)+'].dat->elem_size;')

    code('')
    comm('set up initial pointers and exchange halos if necessary')
    code('int d_m[OPS_MAX_DIM];')
    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_dat':
        code('#ifdef OPS_MPI')
        code('for (int d = 0; d < dim; d++) d_m[d] = args['+str(n)+'].dat->d_m[d] + OPS_sub_dat_list[args['+str(n)+'].dat->index]->d_im[d];')
        code('#else //OPS_MPI')
        code('for (int d = 0; d < dim; d++) d_m[d] = args['+str(n)+'].dat->d_m[d];')
        code('#endif //OPS_MPI')
        code('int base'+str(n)+' = dat'+str(n)+' * 1 * ')
        if prolong[n] == 1:
          code('  ((start[0]/stride_'+str(n)+'[0]) * args['+str(n)+'].stencil->stride[0] - args['+str(n)+'].dat->base[0] - d_m[0]);')
        elif restrict[n] == 1:
          code('  ((start[0]*stride_'+str(n)+'[0]) * args['+str(n)+'].stencil->stride[0] - args['+str(n)+'].dat->base[0] - d_m[0]);')
        else:
          code('  (start[0] * args['+str(n)+'].stencil->stride[0] - args['+str(n)+'].dat->base[0] - d_m[0]);')
        for d in range (1, NDIM):
          line = 'base'+str(n)+' = base'+str(n)+'+ dat'+str(n)+' *\n'
          for d2 in range (0,d):
            line = line + depth*' '+'  args['+str(n)+'].dat->size['+str(d2)+'] *\n'
          code(line[:-1])
          if prolong[n] == 1:
            code('  ((start['+str(d)+']/stride_'+str(n)+'['+str(d)+']) * args['+str(n)+'].stencil->stride['+str(d)+'] - args['+str(n)+'].dat->base['+str(d)+'] - d_m['+str(d)+']);')
          elif restrict[n] == 1:
            code('  ((start['+str(d)+']*stride_'+str(n)+'['+str(d)+']) * args['+str(n)+'].stencil->stride['+str(d)+'] - args['+str(n)+'].dat->base['+str(d)+'] - d_m['+str(d)+']);')
          else:
            code('  (start['+str(d)+'] * args['+str(n)+'].stencil->stride['+str(d)+'] - args['+str(n)+'].dat->base['+str(d)+'] - d_m['+str(d)+']);')
        code('p_a['+str(n)+'] = (char *)args['+str(n)+'].data + base'+str(n)+';')


      elif arg_typ[n] == 'ops_arg_gbl':
        if accs[n] == OPS_READ:
          code('p_a['+str(n)+'] = args['+str(n)+'].data;')
        else:
          code('')
          code('#ifdef OPS_MPI')
          code('p_a['+str(n)+'] = ((ops_reduction)args['+str(n)+'].data)->data + ((ops_reduction)args['+str(n)+'].data)->size * block->index;')
          code('#else //OPS_MPI')
          code('p_a['+str(n)+'] = ((ops_reduction)args['+str(n)+'].data)->data;')
          code('#endif //OPS_MPI')
          code('')
      elif arg_typ[n] == 'ops_arg_idx':
        code('p_a['+str(n)+'] = (char *)arg_idx;')
        code('')
      code('')
    code('')

    code('ops_H_D_exchanges_host(args, '+str(nargs)+');')
    code('ops_halo_exchanges(args,'+str(nargs)+',range);')
    code('ops_H_D_exchanges_host(args, '+str(nargs)+');')
    code('')
    code('ops_timers_core(&c1,&t1);')
    code('OPS_kernels['+str(nk)+'].mpi_time += t1-t2;')
    code('')

    
    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_dat':
        code('xdim'+str(n)+' = args['+str(n)+'].dat->size[0]*args['+str(n)+'].dat->dim;')
        if NDIM==3:
          code('ydim'+str(n)+' = args['+str(n)+'].dat->size[1];')
    code('')

    code('int n_x;')
    
    if not(MULTI_GRID) :    
    ###################### NON-MULTIGRID LOOP EXECUTION ########################
      if NDIM==3:
        FOR('n_z','start[2]','end[2]')
  
      FOR('n_y','start[1]','end[1]')
      code('#pragma novector')
      code('for( n_x=start[0]; n_x<start[0]+((end[0]-start[0])/SIMD_VEC)*SIMD_VEC; n_x+=SIMD_VEC ) {')
      depth = depth+2
  
      comm('call kernel function, passing in pointers to data -vectorised')
      if reduction == 0 and arg_idx == 0:
        code('#pragma simd')
      FOR('i','0','SIMD_VEC')
      text = name+'( '
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
          text = text +' ('+typs[n]+' *)p_a['+str(n)+']+ i*'+str(stride[NDIM*n])
        else:
          text = text +' ('+typs[n]+' *)p_a['+str(n)+']'
        if nargs <> 1 and n != nargs-1:
          text = text + ','
        else:
          text = text +' );\n'
        if n%n_per_line == 2 and n <> nargs-1:
          text = text +'\n          '
      code(text);
      if arg_idx:
        code('arg_idx[0]++;')
      ENDFOR()
      code('')
  
      comm('shift pointers to data x direction')
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_0)*SIMD_VEC;')
  
      ENDFOR()
      code('')
  
      FOR('n_x','start[0]+((end[0]-start[0])/SIMD_VEC)*SIMD_VEC','end[0]')
      comm('call kernel function, passing in pointers to data - remainder')
      text = name+'( '
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
          text = text +' ('+typs[n]+' *)p_a['+str(n)+']'
        else:
          text = text +' ('+typs[n]+' *)p_a['+str(n)+']'
        if nargs <> 1 and n != nargs-1:
          text = text + ','
        else:
          text = text +' );\n'
        if n%n_per_line == 2 and n <> nargs-1:
          text = text +'\n          '
      code(text);
      code('')
  
  
      comm('shift pointers to data x direction')
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_0);')
  
      if arg_idx:
        code('arg_idx[0]++;')
      ENDFOR()
      code('')
  
  
      comm('shift pointers to data y direction')
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
          code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_1);')
      if arg_idx:
        code('#ifdef OPS_MPI')
        for n in range (0,1):
          code('arg_idx['+str(n)+'] = arg_idx_'+str(n)+';')
        code('#else //OPS_MPI')
        for n in range (0,1):
          code('arg_idx['+str(n)+'] = start['+str(n)+'];')
        code('#endif //OPS_MPI')
        code('arg_idx[1]++;')
      ENDFOR()
  
      if NDIM==3:
        comm('shift pointers to data z direction')
        for n in range (0, nargs):
          if arg_typ[n] == 'ops_arg_dat':
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_2);')
  
        if arg_idx:
          code('#ifdef OPS_MPI')
          for n in range (0,2):
            code('arg_idx['+str(n)+'] = sb->decomp_disp['+str(n)+']+start['+str(n)+'];')
          code('#else //OPS_MPI')
          for n in range (0,2):
            code('arg_idx['+str(n)+'] = start['+str(n)+'];')
          code('#endif //OPS_MPI')
          code('arg_idx[2]++;')
        ENDFOR()
    else:
    ######################### MULTIGRID LOOP EXECUTION #########################
      if NDIM==3:
        FOR('n_z','start[2]','end[2]')
  
      FOR('n_y','start[1]','end[1]')
      code('#pragma novector')
      code('for( n_x=start[0]; n_x<end[0]; n_x++ ) {')
      depth = depth+2
  
      comm('call kernel function, passing in pointers to data')
      text = name+'( '
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
          text = text +' ('+typs[n]+' *)p_a['+str(n)+']'
        else:
          text = text +' ('+typs[n]+' *)p_a['+str(n)+']'
        if nargs <> 1 and n != nargs-1:
          text = text + ','
        else:
          text = text +' );\n'
        if n%n_per_line == 2 and n <> nargs-1:
          text = text +'\n          '
      code(text);
      code('')  
  
      comm('shift pointers to data x direction')
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
          if restrict[n] == 1:
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_0) * stride_'+str(n)+'[0];')
          elif prolong[n] == 1:
            #code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_0) * (((global_idx[0]+1) % stride_'+str(n)+'[0] == start[0]% stride_'+str(n)+'[0])?1:0);')
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_0) * (((global_idx[0]+1) % stride_'+str(n)+'[0] == 0)?1:0);')
          else:
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_0);')            
  
      if arg_idx:
        code('arg_idx[0]++;')
      if MULTI_GRID:
        code('global_idx[0]++;')
        
      ENDFOR()
      code('')  
  
      comm('shift pointers to data y direction')
      for n in range (0, nargs):
        if arg_typ[n] == 'ops_arg_dat':
          if restrict[n] == 1:
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_1) * stride_'+str(n)+'[1];')
          elif prolong[n] == 1:
            #IF('(global_idx[1]+1) % stride_'+str(n)+'[1] == start[1] % stride_'+str(n)+'[1]')
            IF('(global_idx[1]+1) % stride_'+str(n)+'[1] == 0')
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_1);')
            ENDIF()
            ELSE()
            code('p_a['+str(n)+']= p_a['+str(n)+'] - (dat'+str(n)+' * off'+str(n)+'_0) * (end_'+str(n)+'[0]-start_'+str(n)+'[0]);')
            ENDIF()
          else:
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_1);')
          
      if arg_idx:
        code('')
        code('#ifdef OPS_MPI')
        for n in range (0,1):
          code('arg_idx['+str(n)+'] = sb->decomp_disp['+str(n)+']+start['+str(n)+'];')
        code('#else //OPS_MPI')
        for n in range (0,1):
          code('arg_idx['+str(n)+'] = start['+str(n)+'];')
        code('#endif //OPS_MPI')
        code('arg_idx[1]++;')
        
      if MULTI_GRID:
        code('')
        code('#ifdef OPS_MPI')
        for n in range (0,1):
          code('global_idx['+str(n)+'] = sb->decomp_disp['+str(n)+']+start['+str(n)+'];')
        code('#else //OPS_MPI')
        for n in range (0,1):
          code('global_idx['+str(n)+'] = start['+str(n)+'];')
        code('#endif //OPS_MPI')
        code('global_idx[1]++;')
        
      ENDFOR()
  
      if NDIM==3:
        comm('shift pointers to data z direction')
        for n in range (0, nargs):
          if arg_typ[n] == 'ops_arg_dat':
            code('p_a['+str(n)+']= p_a['+str(n)+'] + (dat'+str(n)+' * off'+str(n)+'_2);')
  
        if arg_idx:
          code('#ifdef OPS_MPI')
          for n in range (0,2):
            code('arg_idx['+str(n)+'] = sb->decomp_disp['+str(n)+']+start['+str(n)+'];')
          code('#else //OPS_MPI')
          for n in range (0,2):
            code('arg_idx['+str(n)+'] = start['+str(n)+'];')
          code('#endif //OPS_MPI')
          code('arg_idx[2]++;')
        ENDFOR()     
    
    
    code('ops_timers_core(&c2,&t2);')
    code('OPS_kernels['+str(nk)+'].time += t2-t1;')

    code('ops_set_dirtybit_host(args, '+str(nargs)+');')
    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_dat' and (accs[n] == OPS_WRITE or accs[n] == OPS_RW or accs[n] == OPS_INC):
        code('ops_set_halo_dirtybit3(&args['+str(n)+'],range);')

    code('')
    comm('Update kernel record')
    for n in range (0, nargs):
      if arg_typ[n] == 'ops_arg_dat':
        code('OPS_kernels['+str(nk)+'].transfer += ops_compute_transfer(dim, range, &arg'+str(n)+');')
    depth = depth - 2
    code('}')

##########################################################################
#  output individual kernel file
##########################################################################
    if not os.path.exists('./MPI'):
      os.makedirs('./MPI')
    fid = open('./MPI/'+name+'_seq_kernel.cpp','w')
    date = datetime.datetime.now()
    fid.write('//\n// auto-generated by ops.py\n//\n')
    fid.write(file_text)
    fid.close()

# end of main kernel call loop

##########################################################################
#  output one master kernel file
##########################################################################
  depth = 0
  file_text =''
  comm('header')
  if NDIM==3:
    code('#define OPS_3D')
  code('#include "ops_lib_cpp.h"')
  code('#ifdef OPS_MPI')
  code('#include "ops_mpi_core.h"')
  code('#endif')
  if os.path.exists('./user_types.h'):
    code('#include "user_types.h"')
  code('')

  comm(' global constants')
  for nc in range (0,len(consts)):
    if consts[nc]['dim'].isdigit() and int(consts[nc]['dim'])==1:
      code('extern '+consts[nc]['type']+' '+(str(consts[nc]['name']).replace('"','')).strip()+';')
#      code('#pragma acc declare create('+(str(consts[nc]['name']).replace('"','')).strip()+')')
    else:
      if consts[nc]['dim'].isdigit() and consts[nc]['dim'] > 0:
        num = str(consts[nc]['dim'])
        code('extern '+consts[nc]['type']+' '+(str(consts[nc]['name']).replace('"','')).strip()+'['+num+'];')
#        code('#pragma acc declare create('+(str(consts[nc]['name']).replace('"','')).strip()+')')
      else:
        code('extern '+consts[nc]['type']+' *'+(str(consts[nc]['name']).replace('"','')).strip()+';')
#        code('#pragma acc declare create('+(str(consts[nc]['name']).replace('"','')).strip()+')')

  #constants for macros -- now included in teh backend so no need to generate here
  #for i in range(0,20):
  #  code('int xdim'+str(i)+';')
  #code('')

  comm('user kernel files')

  kernel_name_list = []

  for nk in range(0,len(kernels)):
    if kernels[nk]['name'] not in kernel_name_list :
      code('#include "'+kernels[nk]['name']+'_seq_kernel.cpp"')
      kernel_name_list.append(kernels[nk]['name'])

  master = master.split('.')[0]
  fid = open('./MPI/'+master.split('.')[0]+'_seq_kernels.cpp','w')
  fid.write('//\n// auto-generated by ops.py//\n\n')
  fid.write(file_text)
  fid.close()
