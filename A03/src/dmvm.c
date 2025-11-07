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

double dmvm(double *restrict y,
            const double *restrict a,
            double *restrict x,
            int N,
            int Nlocal,
            int iter)
{
  double ts, te;
  int rank, size;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int num = N / size;
  int rest = N % size;

  int upperNeighbor = upper_nbr(rank, size);
  int lowerNeighbor = lower_nbr(rank, size);
  int cs = rows_start_of_rank(rank, size, N);

  ts = MPI_Wtime();
  for (int j = 0; j < iter; j++)
  {
    int currentN = Nlocal;
    int rankCurrent = rank;

    // loop over RHS ring shifts
    for (int rot = 0; rot < size; rot++)
    {

      // local DMVM
      for (int r = 0; r < Nlocal; r++)
      {
        for (int c = cs; c < cs + currentN; c++)
        {
          y[r] = y[r] + a[r * N + c] * x[c - cs];
        }
      }

      // ringshift communication
      cs += currentN;
      if (cs >= N)
      {
        cs = 0; // wrap around
      }

      rankCurrent++;
      if (rankCurrent == size)
      {
        rankCurrent = 0;
      }
      currentN = rows_in_rank(rankCurrent, size, N);

      if (rot != (size - 1))
      {
        // We send upwards towards lower ranks
        MPI_Send(x, num + (rest ? 1 : 0), MPI_DOUBLE, upperNeighbor, 0, MPI_COMM_WORLD);

        MPI_Recv(x, num + (rest ? 1 : 0), MPI_DOUBLE, lowerNeighbor, 0, MPI_COMM_WORLD, &status);

        // MPI_Sendrecv_replace(x,
        //     num + (rest ? 1 : 0),
        //     MPI_DOUBLE,
        //     upperNeighbor,
        //     0,
        //     lowerNeighbor,
        //     0,
        //     MPI_COMM_WORLD,
        //     &status);
      }
    }
  }
  te = MPI_Wtime();

  return te - ts;
}