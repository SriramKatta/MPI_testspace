/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "comm.h"
#include "vtkWriter.h"
#define G(v, i, j, k) v[(k)*imax * jmax + (j)*imax + (i)]
#define L(v, i, j, k) v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

#if defined(_MPI)
// reset fileview for output of string headers
static void resetFileview(VtkOptions* o)
{
    MPI_Offset disp;
    MPI_File_sync(o->fh);
    MPI_Barrier(o->comm.comm);
    MPI_File_get_size(o->fh, &disp);
    MPI_File_set_view(o->fh, disp, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
}
#endif

static double floatSwap(double f)
{
    union {
        double f;
        char b[8];
    } dat1, dat2;

    dat1.f    = f;
    dat2.b[0] = dat1.b[7];
    dat2.b[1] = dat1.b[6];
    dat2.b[2] = dat1.b[5];
    dat2.b[3] = dat1.b[4];
    dat2.b[4] = dat1.b[3];
    dat2.b[5] = dat1.b[2];
    dat2.b[6] = dat1.b[1];
    dat2.b[7] = dat1.b[0];
    return dat2.f;
}

static void writeHeader(VtkOptions* o)
{
#if defined(_MPI)
    // Only rank 0 should write this header to the file.
    if (commIsMaster(&o->comm)) {
        char header[512];
        int len = 0;
        
        len += snprintf(header + len, sizeof(header) - len, "# vtk DataFile Version 3.0\n");
        len += snprintf(header + len, sizeof(header) - len, "PAMPI cfd solver output\n");
        if (o->fmt == ASCII) {
            len += snprintf(header + len, sizeof(header) - len, "ASCII\n");
        } else if (o->fmt == BINARY) {
            len += snprintf(header + len, sizeof(header) - len, "BINARY\n");
        }
        len += snprintf(header + len, sizeof(header) - len, "DATASET STRUCTURED_POINTS\n");
        len += snprintf(header + len, sizeof(header) - len, "DIMENSIONS %d %d %d\n", 
                        o->grid.imax, o->grid.jmax, o->grid.kmax);
        len += snprintf(header + len, sizeof(header) - len, "ORIGIN %f %f %f\n",
                        o->grid.dx * 0.5, o->grid.dy * 0.5, o->grid.dz * 0.5);
        len += snprintf(header + len, sizeof(header) - len, "SPACING %f %f %f\n", 
                        o->grid.dx, o->grid.dy, o->grid.dz);
        len += snprintf(header + len, sizeof(header) - len, "POINT_DATA %d\n", 
                        o->grid.imax * o->grid.jmax * o->grid.kmax);
        
        MPI_File_write(o->fh, header, len, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    resetFileview(o);
#else
    fprintf(o->fh, "# vtk DataFile Version 3.0\n");
    fprintf(o->fh, "PAMPI cfd solver output\n");
    if (o->fmt == ASCII) {
        fprintf(o->fh, "ASCII\n");
    } else if (o->fmt == BINARY) {
        fprintf(o->fh, "BINARY\n");
    }

    fprintf(o->fh, "DATASET STRUCTURED_POINTS\n");
    fprintf(o->fh, "DIMENSIONS %d %d %d\n", o->grid.imax, o->grid.jmax, o->grid.kmax);
    fprintf(o->fh,
        "ORIGIN %f %f %f\n",
        o->grid.dx * 0.5,
        o->grid.dy * 0.5,
        o->grid.dz * 0.5);
    fprintf(o->fh, "SPACING %f %f %f\n", o->grid.dx, o->grid.dy, o->grid.dz);
    fprintf(o->fh, "POINT_DATA %d\n", o->grid.imax * o->grid.jmax * o->grid.kmax);
#endif
}

void vtkOpen(VtkOptions* o, char* problem)
{
    char filename[50];
    snprintf(filename, 50, "%s.vtk", problem);

#if defined(_MPI)
    // Open the vtk file using MPI IO
    MPI_File_open(o->comm.comm, filename, 
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, 
                  MPI_INFO_NULL, &o->fh);
    writeHeader(o);

    if (commIsMaster(&o->comm)) {
        printf("Writing VTK output for %s\n", problem);
    }
#else
    o->fh = fopen(filename, "w");
    writeHeader(o);

    printf("Writing VTK output for %s\n", problem);
#endif
}

static void writeScalar(VtkOptions* o, double* s)
{
    int imax = o->grid.imax;
    int jmax = o->grid.jmax;
    int kmax = o->grid.kmax;

    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                if (o->fmt == ASCII) {
                    fprintf(o->fh, "%f\n", G(s, i, j, k));
                } else if (o->fmt == BINARY) {
                    fwrite((double[1]) { floatSwap(G(s, i, j, k)) },
                        sizeof(double),
                        1,
                        o->fh);
                }
            }
        }
    }
    if (o->fmt == BINARY) fprintf(o->fh, "\n");
}

#if defined(_MPI)
static void commGetOffsets(Comm* c, int offsets[NDIMS], int imaxLocal[], int jmaxLocal[], int kmaxLocal[])
{
    // Gather local sizes from all ranks
    MPI_Allgather(&c->imaxLocal, 1, MPI_INT, imaxLocal, 1, MPI_INT, c->comm);
    MPI_Allgather(&c->jmaxLocal, 1, MPI_INT, jmaxLocal, 1, MPI_INT, c->comm);
    MPI_Allgather(&c->kmaxLocal, 1, MPI_INT, kmaxLocal, 1, MPI_INT, c->comm);
    
    // Calculate offsets based on coordinates in the Cartesian topology
    int coords[NCORDS];
    MPI_Cart_coords(c->comm, c->rank, NDIMS, coords);
    
    // Sum up sizes of ranks before this one in each dimension
    offsets[IDIM] = 0;
    offsets[JDIM] = 0;
    offsets[KDIM] = 0;
    
    for (int r = 0; r < c->size; r++) {
        int rcoords[NCORDS];
        MPI_Cart_coords(c->comm, r, NDIMS, rcoords);
        
        // Add to i offset if same j,k but lower i
        if (rcoords[JCORD] == coords[JCORD] && rcoords[KCORD] == coords[KCORD] && rcoords[ICORD] < coords[ICORD]) {
            offsets[IDIM] += imaxLocal[r];
        }
        // Add to j offset if same i,k but lower j
        if (rcoords[ICORD] == coords[ICORD] && rcoords[KCORD] == coords[KCORD] && rcoords[JCORD] < coords[JCORD]) {
            offsets[JDIM] += jmaxLocal[r];
        }
        // Add to k offset if same i,j but lower k
        if (rcoords[ICORD] == coords[ICORD] && rcoords[JCORD] == coords[JCORD] && rcoords[KCORD] < coords[KCORD]) {
            offsets[KDIM] += kmaxLocal[r];
        }
    }
}
#endif

void vtkScalar(VtkOptions* o, char* name, double* s)
{
#if defined(_MPI)
    if (commIsMaster(&o->comm)) {
        printf("Register scalar %s\n", name);
    }
    
    // 1. Reset fileview before starting MPI IO
    resetFileview(o);
    
    // Write header (only rank 0)
    if (commIsMaster(&o->comm)) {
        char header[128];
        int len = snprintf(header, sizeof(header), "SCALARS %s double 1\nLOOKUP_TABLE default\n", name);
        MPI_File_write(o->fh, header, len, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    // 2. Get offsets from all ranks
    int imaxLocal[o->comm.size];
    int jmaxLocal[o->comm.size];
    int kmaxLocal[o->comm.size];
    int offsets[NDIMS];
    commGetOffsets(&o->comm, offsets, imaxLocal, jmaxLocal, kmaxLocal);
    
    // 3. Sync and get file size for displacement
    MPI_Offset disp;
    MPI_File_sync(o->fh);
    MPI_Barrier(o->comm.comm);
    MPI_File_get_size(o->fh, &disp);
    
    // 4. Create subarray type for global file view
    int globalSizes[NDIMS] = { o->grid.kmax, o->grid.jmax, o->grid.imax };
    int localSizes[NDIMS] = { o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal };
    int starts[NDIMS] = { offsets[KDIM], offsets[JDIM], offsets[IDIM] };
    
    MPI_Datatype fileType;
    MPI_Type_create_subarray(NDIMS, globalSizes, localSizes, starts,
                             MPI_ORDER_C, MPI_DOUBLE, &fileType);
    MPI_Type_commit(&fileType);
    
    // 5. Set file view for this rank
    MPI_File_set_view(o->fh, disp, MPI_DOUBLE, fileType, "native", MPI_INFO_NULL);
    
    // 6. Create subarray for local memory (excluding ghost cells)
    int imaxLocalWithGhost = o->comm.imaxLocal + 2;
    int jmaxLocalWithGhost = o->comm.jmaxLocal + 2;
    int kmaxLocalWithGhost = o->comm.kmaxLocal + 2;
    
    int memorySizes[NDIMS] = { kmaxLocalWithGhost, jmaxLocalWithGhost, imaxLocalWithGhost };
    int memoryStarts[NDIMS] = { 1, 1, 1 };  // Skip ghost cells
    
    MPI_Datatype memType;
    MPI_Type_create_subarray(NDIMS, memorySizes, localSizes, memoryStarts,
                             MPI_ORDER_C, MPI_DOUBLE, &memType);
    MPI_Type_commit(&memType);
    
    // 7. Write data using collective IO
    MPI_File_write_all(o->fh, s, 1, memType, MPI_STATUS_IGNORE);
    
    MPI_Type_free(&fileType);
    MPI_Type_free(&memType);
    
    // Reset fileview and write newline for binary format
    resetFileview(o);
    if (o->fmt == BINARY && commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
    }
#else
    printf("Register scalar %s\n", name);
    if (!isInitialized(o->fh)) return;
    
    fprintf(o->fh, "SCALARS %s double 1\n", name);
    fprintf(o->fh, "LOOKUP_TABLE default\n");
    writeScalar(o, s);
#endif
}

static void writeVector(VtkOptions* o, VtkVector vec)
{
    int imax = o->grid.imax;
    int jmax = o->grid.jmax;
    int kmax = o->grid.kmax;

    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                if (o->fmt == ASCII) {
                    fprintf(o->fh,
                        "%f %f %f\n",
                        G(vec.u, i, j, k),
                        G(vec.v, i, j, k),
                        G(vec.w, i, j, k));
                } else if (o->fmt == BINARY) {
                    fwrite((double[3]) { floatSwap(G(vec.u, i, j, k)),
                               floatSwap(G(vec.v, i, j, k)),
                               floatSwap(G(vec.w, i, j, k)) },
                        sizeof(double),
                        3,
                        o->fh);
                }
            }
        }
    }
    if (o->fmt == BINARY) fprintf(o->fh, "\n");
}

void vtkVector(VtkOptions* o, char* name, VtkVector vec)
{
#if defined(_MPI)
    if (commIsMaster(&o->comm)) {
        printf("Register vector %s\n", name);
    }
    
    // 1. Reset fileview before starting MPI IO
    resetFileview(o);
    
    // Write header (only rank 0)
    if (commIsMaster(&o->comm)) {
        char header[128];
        int len = snprintf(header, sizeof(header), "VECTORS %s double\n", name);
        MPI_File_write(o->fh, header, len, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    // 2. Get offsets from all ranks
    int imaxLocal[o->comm.size];
    int jmaxLocal[o->comm.size];
    int kmaxLocal[o->comm.size];
    int offsets[NDIMS];
    commGetOffsets(&o->comm, offsets, imaxLocal, jmaxLocal, kmaxLocal);
    
    // 3. Sync and get file size for displacement
    MPI_Offset disp;
    MPI_File_sync(o->fh);
    MPI_Barrier(o->comm.comm);
    MPI_File_get_size(o->fh, &disp);
    
    // Create a contiguous type for 3 doubles (u, v, w)
    MPI_Datatype vec3Type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &vec3Type);
    MPI_Type_commit(&vec3Type);
    
    // 4. Create subarray type for global file view (using vec3 as base type)
    int globalSizes[NDIMS] = { o->grid.kmax, o->grid.jmax, o->grid.imax };
    int localSizes[NDIMS] = { o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal };
    int starts[NDIMS] = { offsets[KDIM], offsets[JDIM], offsets[IDIM] };
    
    MPI_Datatype fileType;
    MPI_Type_create_subarray(NDIMS, globalSizes, localSizes, starts,
                             MPI_ORDER_C, vec3Type, &fileType);
    MPI_Type_commit(&fileType);
    
    // 5. Set file view for this rank
    MPI_File_set_view(o->fh, disp, vec3Type, fileType, "native", MPI_INFO_NULL);
    
    // 6. Pack local data (excluding ghost cells) into temporary buffer
    int imaxLocalVal = o->comm.imaxLocal;
    int jmaxLocalVal = o->comm.jmaxLocal;
    int kmaxLocalVal = o->comm.kmaxLocal;
    int imaxLocalWithGhost = imaxLocalVal + 2;
    int jmaxLocalWithGhost = jmaxLocalVal + 2;
    int kmaxLocalWithGhost = kmaxLocalVal + 2;
    
    // Allocate temporary buffer for interleaved vector data
    size_t bufSize = imaxLocalVal * jmaxLocalVal * kmaxLocalVal * 3;
    double* tmpBuf = (double*)malloc(bufSize * sizeof(double));
    
    int idx = 0;
    for (int k = 1; k <= kmaxLocalVal; k++) {
        for (int j = 1; j <= jmaxLocalVal; j++) {
            for (int i = 1; i <= imaxLocalVal; i++) {
                tmpBuf[idx++] = L(vec.u, i, j, k);
                tmpBuf[idx++] = L(vec.v, i, j, k);
                tmpBuf[idx++] = L(vec.w, i, j, k);
            }
        }
    }
    
    // 7. Write data using collective IO
    int count = imaxLocalVal * jmaxLocalVal * kmaxLocalVal;
    MPI_File_write_all(o->fh, tmpBuf, count, vec3Type, MPI_STATUS_IGNORE);
    
    free(tmpBuf);
    MPI_Type_free(&fileType);
    MPI_Type_free(&vec3Type);
    
    // Reset fileview and write newline for binary format
    resetFileview(o);
    if (o->fmt == BINARY && commIsMaster(&o->comm)) {
        MPI_File_write(o->fh, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
    }
#else
    printf("Register vector %s\n", name);
    if (!isInitialized(o->fh)) return;
    
    fprintf(o->fh, "VECTORS %s double\n", name);
    writeVector(o, vec);
#endif
}

void vtkClose(VtkOptions* o)
{
#if defined(_MPI)
    MPI_File_close(&o->fh);
#else
    fclose(o->fh);
    o->fh = NULL;
#endif
}
