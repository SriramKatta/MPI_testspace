/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#ifndef __UTIL_H_
#define __UTIL_H_

#include <mpi.h>

#define HLINE "----------------------------------------------------------------------------\n"

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif
#ifndef ABS
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#endif


#ifdef NDEBUG

#define MPI_CALL(MPIFUNCCALL) \
  {                           \
    MPIFUNCCALL;              \
  }

#else

#define MPI_CALL(MPIFUNCCALL)                                                    \
  {                                                                              \
    int mpi_status = MPIFUNCCALL;                                                \
    if (MPI_SUCCESS != mpi_status)                                               \
    {                                                                            \
      char mpi_error_string[MPI_MAX_ERROR_STRING];                               \
      int mpi_error_string_length = 0;                                           \
      MPI_Error_string(mpi_status, mpi_error_string,                             \
                       &mpi_error_string_length);                                \
      if (NULL != mpi_error_string)                                              \
        fprintf(stderr,                                                          \
                "ERROR: MPI call \"%s\" in line %d of file %s failed "           \
                "with %s "                                                       \
                "(%d).\n",                                                       \
                #MPIFUNCCALL, __LINE__, __FILE__, mpi_error_string, mpi_status); \
      else                                                                       \
        fprintf(stderr,                                                          \
                "ERROR: MPI call \"%s\" in line %d of file %s failed "           \
                "with %d.\n",                                                    \
                #MPIFUNCCALL, __LINE__, __FILE__, mpi_status);                   \
      exit(mpi_status);                                                          \
    }                                                                            \
  }
#endif

#endif // __UTIL_H_
