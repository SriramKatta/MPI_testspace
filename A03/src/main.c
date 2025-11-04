/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "allocate.h"
#include "timing.h"
#include "util.h"

extern double dmvm(double *restrict y, const double *restrict a,
                   const double *restrict x, int Nlocal, int N, int iter);

int main(int argc, char **argv) {
  int rank = 0;
  int size = 1;
  MPI_CALL(MPI_Init(&argc, &argv));
  MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));
  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  size_t bytesPerWord = sizeof(double);
  size_t N = 0;
  size_t iter = 1;
  double *a, *x, *y;
  double t0, t1;
  double walltime;

  if (argc > 2) {
    N = atoi(argv[1]);
    iter = atoi(argv[2]);
  } else {
    if (rank == 0)
      printf("Usage: %s <N> <iter>\n", argv[0]);

    MPI_CALL(MPI_Finalize());
    exit(EXIT_SUCCESS);
  }

  int Nlocal = rows_in_rank(rank, size, N);
  int chunkstart = rows_start_of_rank(rank, size, N);

  a = (double *)allocate(ARRAY_ALIGNMENT, Nlocal * N * bytesPerWord);
  x = (double *)allocate(ARRAY_ALIGNMENT, Nlocal * bytesPerWord);
  y = (double *)allocate(ARRAY_ALIGNMENT, Nlocal * bytesPerWord);

  // initialize arrays
  for (int i = 0; i < Nlocal; i++) {
    x[i] = (double)(i + chunkstart);
    y[i] = 0.0;

    for (int j = 0; j < N; j++) {
      a[i * N + j] = (double)(j + i + chunkstart);
    }
  }

  walltime = dmvm(y, a, x, Nlocal, N, iter);

  // double flops = (double)2.0 * N * N * iter;
  // // # iterations, problem size, flop rate, walltime
  // if(rank == 0)
  // printf("%zu %zu %.2f %.2f\n", iter, N, 1.0E-06 * flops / walltime,
  // walltime);

  free(a);
  free(x);
  free(y);

  MPI_CALL(MPI_Finalize());

  return EXIT_SUCCESS;
}
