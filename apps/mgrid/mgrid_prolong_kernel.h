#ifndef MGRID_PROLONG_KERNELS_H
#define MGRID_PROLONG_KERNELS_H

void mgrid_prolong_kernel(const double *coarse, double *fine, int *idx) {
  //printf("idx = %d, prolong(idx) = %d\n",idx[0], idx[0]/2);
  //printf("idx[0] = %d, idx[1] = %d, idx[0]/2 = %d, idx[1]/2 = %d\n",idx[0], idx[1], idx[0]/2, idx[1]/2);
  //printf("idx[0] = %d, idx[1] = %d, coar %lf\n",
  //        idx[0],idx[1], coarse[OPS_ACC0(0,0)] );
  fine[OPS_ACC1(0,0)] = coarse[OPS_ACC0(0,0)];
}



#endif //MGRID_PROLONG_KERNELS_H
