#define OPS_3D
#include "ops_seq.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {

  //*******************************************************************
  // INITIALISE OPS
  //---------------------------------------
  ops_init(argc, argv, 5);
  ops_printf("Hello world from OPS!\n\n");
  //*******************************************************************

  // THIS IS THE MAIN DIFFERENCE BETWEEN THIS AND THE OTHER PROGRAM
  ops_block block = ops_decl_block_hdf5(3, "grid0", "write_data.h5");

  ops_dat single =
      ops_decl_dat_hdf5(block, 1, "double", "single", "write_data.h5");
  ops_dat multi =
      ops_decl_dat_hdf5(block, 2, "double", "multi", "write_data.h5");

  ops_partition("empty_string_that_does_nothing_yet");
  ops_diagnostic_output();

  ops_fetch_block_hdf5_file(block, "read_data.h5");
  ops_fetch_dat_hdf5_file(multi, "read_data.h5");
  ops_fetch_dat_hdf5_file(single, "read_data.h5");

  //*******************************************************************
  // EXIT OPS AND PRINT TIMING INFO
  //---------------------------------------
  ops_timing_output(stdout);
  ops_printf("\nSucessful exit from OPS!\n");
  ops_exit();
  //*******************************************************************
}
