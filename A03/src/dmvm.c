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
               const double *restrict x, int rows, int cols)
{
  for (int r = 0; r < rows; r++)
  {
    for (int c = 0; c < cols; c++)
    {
      y[r] += a[r * cols + c] * x[c];
    }
  }
#ifdef CHECK
  {
    double sum = 0.0;

    for (int i = 0; i < N; i++)
    {
      sum += y[i];
      y[i] = 0.0;
    }
    fprintf(stderr, "Sum: %f\n", sum);
  }
#endif
}

double dmvm(double *restrict y, const double *restrict a,
            const double *restrict x, int rows, int cols, int iter)
{
  double ts, te;
  int size = 0;
  int rank = 1;

  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));

  int Nlocal = rows_in_rank(rank, size, cols);
  int row_start = rows_start_of_rank(rank, size, cols);
  int Ncurrent = Nlocal;
  int lnbr = lower_nbr(rank, size);
  int unbr = upper_nbr(rank, size);

  ts = getTimeStamp();
  for (int j = 0; j < iter; j++)
  {
    for (int rot = 0; rot < size; ++rot)
    {
      for (int r = 0; r < row_start + Ncurrent; r++)
      {
        for (int c = 0; c < row_start + Nlocal; c++)
        {
          y[r] += a[r * cols + c] * x[c - row_start];
        }
      }

      row_start += Nlocal;
      if (row_start >= cols)
      {
        row_start = 0;
      }
      int lrank = lower_nbr(rank, size);

      Ncurrent = rows_in_rank(rank, size, cols);

      // communicate(x, cols, rank, size);
      // funtion to perform ring shift of x vector
      if (rot != size - 1)
      {
        MPI_CALL(MPI_Sendrecv_replace(x, (cols / size) + (cols % size) ? 1 : 0, MPI_DOUBLE, unbr, 0, lnbr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE));
      }
    }
  }
  te = getTimeStamp();

  return te - ts;
}