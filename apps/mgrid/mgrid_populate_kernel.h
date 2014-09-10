#ifndef MGRID_KERNELS_H
#define MGRID_KERNELS_H

void mgrid_populate_kernel(double *val, int *idx) {
  val[OPS_ACC0(0,0)] = (double)(idx[0]+20*idx[1]);
}

void mgrid_restrict_kernel(const double *fine, double *coarse, int *idx) {
  printf("idx[0] %d, idx[1] %d , coarse[OPS_ACC1(0,0)] %lf, fine[OPS_ACCR0(-1,0)] %lf, fine[OPS_ACCR0(0,0)] %lf, fine[OPS_ACCR0(1,0)] %lf\n",
         idx[0], idx[1], coarse[OPS_ACC1(0,0)], fine[OPS_ACCR0(-1,0)], fine[OPS_ACCR0(0,0)], fine[OPS_ACCR0(1,0)]);
}

void mgrid_prolong_kernel(const double *coarse, double *fine, int *idx) {
  printf("idx[0] %d, idx[1] %d , fine[OPS_ACC1(0,0)] %lf, coarse[OPS_ACCR0(-1,0)] %lf, coarse[OPS_ACCR0(0,0)] %lf, coarse[OPS_ACCR0(1,0)] %lf\n",
         idx[0], idx[1], fine[OPS_ACC1(0,0)], coarse[OPS_ACCR0(-1,0)], coarse[OPS_ACCR0(0,0)], coarse[OPS_ACCR0(1,0)]);
}


#endif //MGRID_KERNELS_H
