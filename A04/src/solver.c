/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"

#include "allocate.h"
#include "parameter.h"
#include "solver.h"
#include "util.h"

#define PI 3.14159265358979323846
#define P(i, j) p[(j) * (imax + 2) + (i)]
#define RHS(i, j) rhs[(j) * (imax + 2) + (i)]

static int sizeOfRank(int rank, int size, int N)
{
    int base = N / size;
    int rest = N % size;
    return base + ((rest > rank) ? 1 : 0);
}

static void exchange(Solver *solver)
{
    MPI_Request requests[4] = {MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL};

    double *p = solver->p;
    int imax = solver->imax;

    /* exchange ghost cells with top neighbor */
    if (solver->rank + 1 < solver->size)
    {
        int top = solver->rank + 1;
        double *src = &P(1, solver->jmaxLocal);
        double *dst = &P(1, solver->jmaxLocal + 1);
        //           col ^, ^ row

        MPI_CALL(MPI_Isend(src, solver->imax, MPI_DOUBLE, top, 1, MPI_COMM_WORLD, &requests[0]));
        MPI_CALL(MPI_Irecv(dst, solver->imax, MPI_DOUBLE, top, 2, MPI_COMM_WORLD, &requests[1]));
    }

    /* exchange ghost cells with bottom neighbor */
    if (solver->rank > 0)
    {
        int bottom = solver->rank - 1;
        double *src = &P(1, 1);
        double *dst = &P(1, 0);
        //           col ^, ^ row

        MPI_CALL(MPI_Isend(src, solver->imax, MPI_DOUBLE, bottom, 2, MPI_COMM_WORLD, &requests[2]));
        MPI_CALL(MPI_Irecv(dst, solver->imax, MPI_DOUBLE, bottom, 1, MPI_COMM_WORLD, &requests[3]));
    }

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
}

void getResult(Solver *solver)
{
    double *pfull = NULL;
    int *rcvCounts, *displs;

    if (solver->rank == 0)
    {
        pfull = (double *)allocate(64, (solver->imax + 2) * (solver->jmax + 2) * sizeof(double));
        rcvCounts = (int *)malloc(solver->size * sizeof(int));
        displs = (int *)malloc(solver->size * sizeof(int));
        rcvCounts[0] = solver->jmaxLocal * (solver->imax + 2);
        displs[0] = 0;
        int cursor = rcvCounts[0];

        for (int i = 1; i < solver->size; i++)
        {
            rcvCounts[i] = sizeOfRank(i, solver->size, solver->jmax) * (solver->imax + 2);
            displs[i] = cursor;
            cursor += rcvCounts[i];
        }
    }

    int cnt = solver->jmaxLocal * (solver->imax + 2);
    double *sendbuffer = solver->p + (solver->imax + 2); // to skip the bottom row since its a ghost layer
    MPI_CALL(MPI_Gatherv(sendbuffer,
                         cnt,
                         MPI_DOUBLE,
                         pfull,
                         rcvCounts,
                         displs,
                         MPI_DOUBLE,
                         0,
                         MPI_COMM_WORLD));
    if (solver->rank == 0)
    {
        writeResult(solver, pfull, "p.dat");
        free(pfull);
    }
}

void initSolver(Solver *solver, Parameter *params, int problem)
{
    MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &(solver->rank)));
    MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &(solver->size)));
    solver->imax = params->imax;
    solver->jmax = params->jmax;
    solver->jmaxLocal = sizeOfRank(solver->rank, solver->size, solver->jmax);

    solver->dx = params->xlength / params->imax;
    solver->dy = params->ylength / params->jmax;

    int jstart = 0;
    for (int r = 0; r < solver->rank; ++r)
    {
        jstart += sizeOfRank(r, solver->size, solver->jmax);
    }
    solver->ys = jstart * solver->dy;

    solver->eps = params->eps;
    solver->omega = params->omg;
    solver->itermax = params->itermax;

#ifdef DEBUG
    printf("RANK %d: imaxLocal : %d, jmaxLocal : %d ys: %lf\n",
           solver->rank,
           solver->imax,
           solver->jmaxLocal,
           solver->ys);
#endif

    int imax = solver->imax;
    int jmax = solver->jmax;
    int jmaxlocal = solver->jmaxLocal;
    // adapt for MPI case
    size_t bytesize = (imax + 2) * (jmaxlocal + 2) * sizeof(double);
    solver->p = (double *)allocate(64, bytesize);
    solver->rhs = (double *)allocate(64, bytesize);

    double dx = solver->dx;
    double dy = solver->dy;
    double *p = solver->p;
    double *rhs = solver->rhs;

    // adapted for MPI case
    for (int j = 0; j < jmaxlocal + 2; j++)
    {
        for (int i = 0; i < imax + 2; i++)
        {
            P(i, j) = sin(2.0 * PI * i * dx * 2.0) + sin(2.0 * PI * (j * dy + (solver->ys)) * 2.0);
        }
    }

    if (problem == 2)
    {
        for (int j = 0; j < jmaxlocal + 2; j++)
        {
            for (int i = 0; i < imax + 2; i++)
            {
                RHS(i, j) = sin(2.0 * PI * i * dx);
            }
        }
    }
    else
    {
        for (int j = 0; j < jmaxlocal + 2; j++)
        {
            for (int i = 0; i < imax + 2; i++)
            {
                RHS(i, j) = 0.0;
            }
        }
    }
}

// makes it easier to implement cache blocking to better improve performance
double solver_core(Solver *solver)
{
    double res = 0.0;

    int imax = solver->imax;
    int jmaxlocal = solver->jmaxLocal;
    double *p = solver->p;
    double *rhs = solver->rhs;
    double dx2 = solver->dx * solver->dx;
    double dy2 = solver->dy * solver->dy;
    double idx2 = 1.0 / dx2;
    double idy2 = 1.0 / dy2;
    double factor = solver->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);

    double L2cachebytes = 1.25 * 1000.0 * 1000.0; // cache size in bytes
    int collimit = L2cachebytes / 48;
    collimit = MIN(collimit, imax + 1);

    for (int colstart = 1; colstart < imax + 1; colstart += collimit)
    {
        int colend = MIN(colstart + collimit, imax + 1);
        // adapt for mpi
        for (int j = 1; j < jmaxlocal + 1; j++)
        {
            for (int i = colstart; i < colend; ++i)
            {
                double r = RHS(i, j) -
                           ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                            (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

                P(i, j) -= (factor * r);
                res += (r * r);
            }
        }
    }
    return res;
}

double solver_RB_core(Solver *solver, COLOUR colour)
{
    double res = 0.0;

    int imax = solver->imax;
    int jmaxlocal = solver->jmaxLocal;
    double *p = solver->p;
    double *rhs = solver->rhs;
    double dx2 = solver->dx * solver->dx;
    double dy2 = solver->dy * solver->dy;
    double idx2 = 1.0 / dx2;
    double idy2 = 1.0 / dy2;
    double factor = solver->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);

    int step = ((colour == RED) ? 0 : 1);
    int colstart = step + 1;
    for (int j = 1; j < jmaxlocal + 1; j++)
    {
        for (int i = colstart; i < imax + 1; i += 2)
        {
            double r = RHS(i, j) -
                       ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                        (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

            P(i, j) -= (factor * r);
            res += (r * r);
        }
        colstart = 3 - colstart; // swap between 1 and 2
    }

    return res;
}

void apply_bdy_condition(Solver *solver)
{
    int imax = solver->imax;
    double *p = solver->p;
    int jmaxlocal = solver->jmaxLocal;
    // adapt for mpi // boundary conditions
    for (int i = 1; i < imax + 1; i++)
    {
        if (solver->rank == 0)
        {
            P(i, 0) = P(i, 1);
        }
        if (solver->rank == (solver->size - 1))
        {
            P(i, jmaxlocal + 1) = P(i, jmaxlocal);
        }
    }

    for (int j = 1; j < jmaxlocal + 1; j++)
    {
        P(0, j) = P(1, j);
        P(imax + 1, j) = P(imax, j);
    }
}

void solve(Solver *solver)
{
    int imax = solver->imax;
    int jmax = solver->jmax;
    int jmaxlocal = solver->jmaxLocal;
    double eps = solver->eps;
    int itermax = solver->itermax;
    double *p = solver->p;
    double epssq = eps * eps;
    int it = 0;
    double res = eps + 1.0;

    MPI_Request req;

    while ((res >= epssq) && (it < itermax))
    {
        exchange(solver);
#ifdef SOR_RB_SOLVER
        res += solver_RB_core(solver, RED);
        exchange(solver);
        res += solver_RB_core(solver, BLACK);
#else
        res = solver_core(solver);
#endif

        MPI_CALL(MPI_Iallreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req));

        apply_bdy_condition(solver);

        // termination condition collection
        MPI_CALL(MPI_Wait(&req, MPI_STATUS_IGNORE));
        res = res / (double)(imax * jmax);
#ifdef DEBUG
        if (solver->rank == 0)
        {
            printf("%d Residuum: %e\n", it, res);
        }
#endif
        it++;
    }

    if (solver->rank == 0)
    {
        printf("Solver took %d iterations to reach %f using omega=%f\n",
               it,
               sqrt(res),
               solver->omega);
    }
}


void writeResult(Solver *solver, double *p, char *filename)
{
    int imax = solver->imax;
    int jmax = solver->jmax;

    FILE *fp;
    fp = fopen(filename, "w");

    if (fp == NULL)
    {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < jmax + 2; j++)
    {
        for (int i = 0; i < imax + 2; i++)
        {
            fprintf(fp, "%f ", P(i, j));
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}
