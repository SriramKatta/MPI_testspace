/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>

#define PI 3.14159265

#include "timing.h"

#define DX 1E-9;

#define MPI_CALL(call)                                                         \
  {                                                                            \
    int mpi_status = call;                                                     \
    if (MPI_SUCCESS != mpi_status) {                                           \
      char mpi_error_string[MPI_MAX_ERROR_STRING];                             \
      int mpi_error_string_length = 0;                                         \
      MPI_Error_string(mpi_status, mpi_error_string,                           \
                       &mpi_error_string_length);                              \
      if (NULL != mpi_error_string)                                            \
        fprintf(stderr,                                                        \
                "ERROR: MPI call \"%s\" in line %d of file %s failed "         \
                "with %s "                                                     \
                "(%d).\n",                                                     \
                #call, __LINE__, __FILE__, mpi_error_string, mpi_status);      \
      else                                                                     \
        fprintf(stderr,                                                        \
                "ERROR: MPI call \"%s\" in line %d of file %s failed "         \
                "with %d.\n",                                                  \
                #call, __LINE__, __FILE__, mpi_status);                        \
      exit(mpi_status);                                                        \
    }                                                                          \
  }

double integrate(double, double);

int main(int argc, char **argv) {
  double wcs, wce;
  double Pi = 0.0;

  int rank = 0;
  int size = 1;
  MPI_CALL(MPI_Init(&argc, &argv));
  MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));
  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  double chunk = 1.0 / size;
  double start = rank * chunk;
  double end = start + chunk;

  
  wcs = getTimeStamp();
  double Pi_local = integrate(start, end);
  // printf("rank %d | start %3.4lf | end %3.4lf | Pi local %3.4lf\n", rank, start, end, Pi_local);
  
  // MPI_CALL(MPI_Barrier(MPI_COMM_WORLD));

  if (rank == 0) {
    Pi = Pi_local;
    for (size_t i = 1; i < size; i++) {
      double recv_pi = 0.0;
      MPI_CALL(MPI_Recv(&recv_pi, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE));
      Pi += recv_pi;
    }
  } else {
    MPI_CALL(MPI_Send(&Pi_local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD));
  }
  Pi *= 4.0;
  wce = getTimeStamp();
  if (rank == 0)
    printf("Pi=%1.8lf in %.3lf s \n", Pi, wce - wcs);

  MPI_CALL(MPI_Finalize());
  return EXIT_SUCCESS;
}

double integrate(double a, double b) {

  /*

          Your logic to integrate between given interval a to b.
          Declare SLICES here and calculate delta x using a, b and SLICES.
          Iterate over number of slices, calculate the area and sum them.
          Return sum * delta x to get the value of PI.

  */
  double res = 0.0;
  while (a < b) {
    res += sqrt(1 - a * a);
    a += DX;
  }

  return res * DX;
}
