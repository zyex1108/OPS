#ifndef MGRID_POPULATE_KERNEL_H
#define MGRID_POPULATE_KERNEL_H

void mgrid_populate_kernel(double *val, int *idx) {
  val[OPS_ACC0(0,0)] = (double)(idx[0]+20*idx[1]);
}


#endif //MGRID_POPULATE_KERNEL_H
