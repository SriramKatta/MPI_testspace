/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
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

#include "likwid-marker.h"
#include "parameter.h"
#include "solver.h"
#include "timing.h"
#include "util.h"

int main(int argc, char **argv)
{

    MPI_CALL(MPI_Init(&argc, &argv));
    int rank, size;
    MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &size));

    Parameter params;
    Solver solver;
    initParameter(&params);

    if (argc < 2)
    {
        if (rank == 0)
        {
            printf("Usage: %s <configFile>\n", argv[0]);
        }
        MPI_CALL(MPI_Finalize());
        exit(EXIT_SUCCESS);
    }

    readParameter(&params, argv[1]);

    if (rank == 0)
    {
        printParameter(&params);
    }

    initSolver(&solver, &params, 2);
    double startTime = getTimeStamp();
    solve(&solver);
    double endTime = getTimeStamp();
    // writeResult(&solver, "p.dat");
    getResult(&solver);
    
    if (rank == 0)
    {
        printf("Walltime %.2fs\n", endTime - startTime);
    }

    MPI_CALL(MPI_Finalize());
    return EXIT_SUCCESS;
}
