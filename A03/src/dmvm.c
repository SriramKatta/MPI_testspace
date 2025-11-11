/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "allocate.h"
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

inline void swap_buffer(double **buff1, double **buff2)
{
  double *tmpbuff = *buff1;
  *buff1 = *buff2;
  *buff2 = tmpbuff;
}

double dmvm(double *restrict y,
            const double *restrict a,
            double *restrict x,
            int N,
            int Nlocal,
            int iter)
{
  double *xbuff[2];
  xbuff[0] = x;
#ifdef NB_COMMUICATION
  MPI_Request requests[2];
  double *xtbuff = (double *)allocate(ARRAY_ALIGNMENT, (Nlocal + 1) * sizeof(double));
  xbuff[1] = xtbuff;
#endif

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

    for (int rot = 0; rot < size; rot++)
    {
#ifdef NB_COMMUICATION
      if (rot != (size - 1))
      {
        MPI_CALL(MPI_Isend(xbuff[0], num + (rest ? 1 : 0), MPI_DOUBLE,
                           upperNeighbor, 0, MPI_COMM_WORLD, &requests[0]));

        MPI_CALL(MPI_Irecv(xbuff[1], num + (rest ? 1 : 0), MPI_DOUBLE,
                           lowerNeighbor, 0, MPI_COMM_WORLD, &requests[1]));
      }

#endif

      dmvm_core(Nlocal, cs, currentN, y, a, N, x);

      // set up comm
      cs += currentN;
      if (cs >= N)
      {
        cs = 0;
      }

      rankCurrent = lower_nbr(rankCurrent, size);

      currentN = rows_in_rank(rankCurrent, size, N);

#ifdef NB_COMMUICATION
      if (rot != (size - 1))
      {
        MPI_Waitall(2, requests, MPI_STATUS_IGNORE);
        swap_buffer(&xbuff[0], &xbuff[1]);
      }
#else
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
#endif
    }
  }
  double te = MPI_Wtime();

  return te - ts;
}
