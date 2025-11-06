/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include "timing.h"
#include "util.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void communicate(double *restrict x, int N, int rank, int size)
{
  int urank = upper_nbr(rank, size);
  int lrank = lower_nbr(rank, size);
  int count_curr_rank = rows_in_rank(rank, size, N);
  int count_upp_rank = rows_in_rank(urank, size, N);
  int count_low_rank = rows_in_rank(lrank, size, N);

  if (rank == 0)
  {
    MPI_CALL(MPI_Send(x, count_curr_rank, MPI_DOUBLE, lrank, 0, MPI_COMM_WORLD))
    MPI_CALL(MPI_Recv(x, count_upp_rank, MPI_DOUBLE, urank, 0, MPI_COMM_WORLD,
                      MPI_STATUS_IGNORE));
  }
  else
  {
    MPI_CALL(MPI_Recv(x, count_upp_rank, MPI_DOUBLE, urank, 0, MPI_COMM_WORLD,
                      MPI_STATUS_IGNORE));
    MPI_CALL(MPI_Send(x, count_curr_rank, MPI_DOUBLE, lrank, 0, MPI_COMM_WORLD))
  }
}

void dmvm_core(double *restrict y, const double *restrict a,
               const double *restrict x, int rows, int cols, int colstart)
{
  for (int r = 0; r < rows; r++)
  {
    for (int c = colstart; c < colstart + cols; c++)
    {
      y[r] += a[r * cols + c] * x[c - colstart];
    }
  }
}

double dmvm(double *restrict y, const double *restrict a,
            const double *restrict x, int N, int iter)
{
  double ts, te;
  int size = 0;
  int rank = 1;

  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));

  int Nlocal = rows_in_rank(rank, size, N);
  int col_start = rows_start_of_rank(rank, size, N);
  int lnbr = lower_nbr(rank, size);
  int unbr = upper_nbr(rank, size);

  ts = MPI_Wtime();
  for (int j = 0; j < iter; j++)
  {
    int Ncurrent = Nlocal;
    int rankcurrent = rank;

    for (int roti = 0; roti < size; roti++)
    {
      dmvm_core(y, a, x, N, Ncurrent, col_start);

      col_start += Ncurrent;
      if (col_start > N)
      {
        col_start = 0;
      }

      rankcurrent = lower_nbr(rankcurrent, size);

      Ncurrent = rows_in_rank(rankcurrent, size, N);

      if (roti != size - 1)
      {
        MPI_CALL(MPI_Sendrecv_replace((void *)x,
                                      (N / size) + ((N % size) ? 1 : 0),
                                      MPI_DOUBLE,
                                      unbr, 0,
                                      lnbr, 0,
                                      MPI_COMM_WORLD,
                                      MPI_STATUS_IGNORE));
      }
    }
  }
  te = MPI_Wtime();

  return te - ts;
}