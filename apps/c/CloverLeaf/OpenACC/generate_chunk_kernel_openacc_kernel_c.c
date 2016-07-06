//
// auto-generated by ops.py
//
#include "./OpenACC/clover_leaf_common.h"

#undef OPS_GPU

int xdim0_generate_chunk_kernel;
int xdim1_generate_chunk_kernel;
int xdim2_generate_chunk_kernel;
int xdim3_generate_chunk_kernel;
int xdim4_generate_chunk_kernel;
int xdim5_generate_chunk_kernel;
int xdim6_generate_chunk_kernel;
int xdim7_generate_chunk_kernel;


#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6
#undef OPS_ACC7


#define OPS_ACC0(x,y) (x+xdim0_generate_chunk_kernel*(y))
#define OPS_ACC1(x,y) (x+xdim1_generate_chunk_kernel*(y))
#define OPS_ACC2(x,y) (x+xdim2_generate_chunk_kernel*(y))
#define OPS_ACC3(x,y) (x+xdim3_generate_chunk_kernel*(y))
#define OPS_ACC4(x,y) (x+xdim4_generate_chunk_kernel*(y))
#define OPS_ACC5(x,y) (x+xdim5_generate_chunk_kernel*(y))
#define OPS_ACC6(x,y) (x+xdim6_generate_chunk_kernel*(y))
#define OPS_ACC7(x,y) (x+xdim7_generate_chunk_kernel*(y))

//user function
inline
void generate_chunk_kernel( const double *vertexx, const double *vertexy,
                     double *energy0, double *density0,
                     double *xvel0,  double *yvel0,
                     const double *cellx, const double *celly) {

  double radius, x_cent, y_cent;
  int is_in = 0;
  int is_in2 = 0;


  energy0[OPS_ACC2(0,0)]= states[0].energy;
  density0[OPS_ACC3(0,0)]= states[0].density;
  xvel0[OPS_ACC4(0,0)]=states[0].xvel;
  yvel0[OPS_ACC5(0,0)]=states[0].yvel;

  for(int i = 1; i<number_of_states; i++) {

    x_cent=states[i].xmin;
    y_cent=states[i].ymin;
    is_in = 0;
    is_in2 = 0;

    if (states[i].geometry == g_rect) {
      for (int i1 = -1; i1 <= 0; i1++) {
        for (int j1 = -1; j1 <= 0; j1++) {
          if(vertexx[OPS_ACC0(1+i1,0)] >= states[i].xmin  && vertexx[OPS_ACC0(0+i1,0)] < states[i].xmax) {
            if(vertexy[OPS_ACC1(0,1+j1)] >= states[i].ymin && vertexy[OPS_ACC1(0,0+j1)] < states[i].ymax) {
              is_in = 1;
            }
          }
        }
      }
      if(vertexx[OPS_ACC0(1,0)] >= states[i].xmin  && vertexx[OPS_ACC0(0,0)] < states[i].xmax) {
        if(vertexy[OPS_ACC1(0,1)] >= states[i].ymin && vertexy[OPS_ACC1(0,0)] < states[i].ymax) {
          is_in2 = 1;
        }
      }
      if (is_in2) {
        energy0[OPS_ACC2(0,0)] = states[i].energy;
        density0[OPS_ACC3(0,0)] = states[i].density;
      }
      if (is_in) {
        xvel0[OPS_ACC4(0,0)] = states[i].xvel;
        yvel0[OPS_ACC5(0,0)] = states[i].yvel;
      }
    }
    else if(states[i].geometry == g_circ) {
      for (int i1 = -1; i1 <= 0; i1++) {
        for (int j1 = -1; j1 <= 0; j1++) {
          radius = sqrt ((cellx[OPS_ACC6(i1,0)] - x_cent) * (cellx[OPS_ACC6(i1,0)] - x_cent) +
                     (celly[OPS_ACC7(0,j1)] - y_cent) * (celly[OPS_ACC7(0,j1)] - y_cent));
          if (radius <= states[i].radius) {
            is_in = 1;
          }
        }
      }
      if (radius <= states[i].radius) is_in2 = 1;

      if (is_in2) {
        energy0[OPS_ACC2(0,0)] = states[i].energy;
        density0[OPS_ACC3(0,0)] = states[i].density;
      }

      if (is_in) {
        xvel0[OPS_ACC4(0,0)] = states[i].xvel;
        yvel0[OPS_ACC5(0,0)] = states[i].yvel;
      }
    }
    else if(states[i].geometry == g_point) {
      for (int i1 = -1; i1 <= 0; i1++) {
        for (int j1 = -1; j1 <= 0; j1++) {
          if(vertexx[OPS_ACC0(i1,0)] == x_cent && vertexy[OPS_ACC1(0,j1)] == y_cent) {
            is_in = 1;
          }
        }
      }
      if(vertexx[OPS_ACC0(0,0)] == x_cent && vertexy[OPS_ACC1(0,0)] == y_cent)
        is_in2 = 1;

      if (is_in2) {
        energy0[OPS_ACC2(0,0)] = states[i].energy;
        density0[OPS_ACC3(0,0)] = states[i].density;
      }

      if (is_in) {
        xvel0[OPS_ACC4(0,0)] = states[i].xvel;
        yvel0[OPS_ACC5(0,0)] = states[i].yvel;
      }
    }
  }
}


#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6
#undef OPS_ACC7



void generate_chunk_kernel_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  double *p_a3,
  double *p_a4,
  double *p_a5,
  double *p_a6,
  double *p_a7,
  int x_size, int y_size) {
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a0,p_a1,p_a2,p_a3,p_a4,p_a5,p_a6,p_a7)
  #pragma acc loop
  #endif
  for ( int n_y=0; n_y<y_size; n_y++ ){
    #ifdef OPS_GPU
    #pragma acc loop
    #endif
    for ( int n_x=0; n_x<x_size; n_x++ ){
      generate_chunk_kernel(  p_a0 + n_x*1*1 + n_y*xdim0_generate_chunk_kernel*0*1,
           p_a1 + n_x*0*1 + n_y*xdim1_generate_chunk_kernel*1*1, p_a2 + n_x*1*1 + n_y*xdim2_generate_chunk_kernel*1*1,
           p_a3 + n_x*1*1 + n_y*xdim3_generate_chunk_kernel*1*1, p_a4 + n_x*1*1 + n_y*xdim4_generate_chunk_kernel*1*1,
           p_a5 + n_x*1*1 + n_y*xdim5_generate_chunk_kernel*1*1, p_a6 + n_x*1*1 + n_y*xdim6_generate_chunk_kernel*0*1,
           p_a7 + n_x*0*1 + n_y*xdim7_generate_chunk_kernel*1*1 );

    }
  }
}
