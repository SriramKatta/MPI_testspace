/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <mpi.h>

#include "allocate.h"
#include "util.h"

extern double dmvm(
    double *restrict y,
    const double *restrict a,
    const double *restrict x,
    int N,
    int Nlocal,
    int iter);

int main(int argc, char **argv)
{
  size_t bytesPerWord = sizeof(double);
  size_t N = 0;
  size_t Nlocal = 0;
  size_t iter = 1;
  int rank, size;

  MPI_CALL(MPI_Init(&argc, &argv));
  MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));
  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  if (argc != 3)
  {
    if (rank == 0)
    {
      printf("Usage: %s <N> <iter>\n", argv[0]);
    }
    MPI_CALL(MPI_Finalize());
    exit(EXIT_SUCCESS);
  }

  N = atoi(argv[1]);
  iter = atoi(argv[2]);
  Nlocal = rows_in_rank(rank, size, N);

  double *a = (double *)allocate(ARRAY_ALIGNMENT, N * Nlocal * bytesPerWord);
  double *x = (double *)allocate(ARRAY_ALIGNMENT, (Nlocal + 1) * bytesPerWord);
  double *y = (double *)allocate(ARRAY_ALIGNMENT, (Nlocal + 1) * bytesPerWord);

  // initialize arrays
  int cs = rows_start_of_rank(rank, size, N);

  for (int i = 0; i < Nlocal; i++)
  {
    x[i] = (double)(cs + i);
    y[i] = 0.0;

    for (int j = 0; j < N; j++)
    {
      a[i * N + j] = (double)(j + cs + i);
    }
  }

  double walltime = dmvm(y, a, x, N, Nlocal, iter);

  double flops = (double)2 * N * N * iter;
  if (rank == 0)
  {
    printf("RES | %zu %d %zu %.2f %.2f\n", iter, size, N, 1.0E-06 * flops / walltime,
           walltime);
  }

  MPI_CALL(MPI_Finalize());
  return EXIT_SUCCESS;
}
