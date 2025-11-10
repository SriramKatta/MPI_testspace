/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

inline void dmvm_core(int Nlocal, int cs, int currentN,
                      double *restrict y, const double *restrict a, int N,
                      double *restrict x)
{
  for (int r = 0; r < Nlocal; r++)
  {
    for (int c = cs; c < cs + currentN; c++)
    {
      y[r] = y[r] + a[r * N + c] * x[c - cs];
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

  int rank, size;

  MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));
  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  int num = N / size;
  int rest = N % size;

  int upperNeighbor = upper_nbr(rank, size);
  int lowerNeighbor = lower_nbr(rank, size);

  int cs = rows_start_of_rank(rank, size, N);

  double ts = MPI_Wtime();
  for (int j = 0; j < iter; j++)
  {
    int currentN = Nlocal;
    int rankCurrent = rank;

    // loop over RHS ring shifts
    for (int rot = 0; rot < size; rot++)
    {

      // local DMVM
      dmvm_core(Nlocal, cs, currentN, y, a, N, x);

      cs += currentN;
      if (cs >= N)
      {
        cs = 0;
      }

      rankCurrent = lower_nbr(rankCurrent, size);

      currentN = rows_in_rank(rankCurrent, size, N);

      if (rot != (size - 1))
      {
        if (rank % 2) // to prevent dead locking
        {
          MPI_CALL(MPI_Send(x, num + (rest ? 1 : 0), MPI_DOUBLE,
                            upperNeighbor, 0, MPI_COMM_WORLD));

          MPI_CALL(MPI_Recv(x, num + (rest ? 1 : 0), MPI_DOUBLE,
                            lowerNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE));
        }
        else
        {
          MPI_CALL(MPI_Recv(x, num + (rest ? 1 : 0), MPI_DOUBLE,
                            lowerNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE));

          MPI_CALL(MPI_Send(x, num + (rest ? 1 : 0), MPI_DOUBLE,
                            upperNeighbor, 0, MPI_COMM_WORLD));
        }
      }
    }
  }
  double te = MPI_Wtime();

  return te - ts;
}
