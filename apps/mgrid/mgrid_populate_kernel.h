#ifndef MGRID_POPULATE_KERNELS_H
#define MGRID_POPULATE_KERNELS_H

void mgrid_populate_kernel(double *val, int *idx) {
  val[OPS_ACC0(0,0)] = (double)(idx[0]+5*idx[1]);
}





#endif //MGRID_POPULATE_KERNELS_H
