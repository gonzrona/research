#include "../headers/structs.h"
#include "../headers/prototypes.h"

void create_columns(double _Complex *send, double _Complex *receive, System sys) {
    MPI mpi = sys.mpi; int i,j,Nx = sys.lat.Nx;
    for (i=0; i<mpi.Int_x[mpi.rank]; i++) {
        for (j=mpi.start_y; j<mpi.end_y; j++) {
            receive[ind(i,j,mpi.Int_x[mpi.rank])] = send[ind(i+mpi.start_x,j-mpi.start_y,Nx)];
        }
    }
    for (i=0; i<mpi.world_size; i++) {
        if (i != mpi.rank) {
            MPI_Request request; MPI_Status status;
            MPI_Isend(&send[start_index(mpi.Int_x,i)], 1, mpi.columns[i], i, 0, MPI_COMM_WORLD,&request);
            MPI_Recv(&receive[start_index(mpi.Int_y,i)*mpi.Int_x[mpi.rank]], 1, mpi.rows[i], i, 0, MPI_COMM_WORLD, &status);
            MPI_Wait(&request, &status);
        }
    }
}

void create_rows(double _Complex *send, double _Complex *receive, System sys) {
    MPI mpi = sys.mpi; int i,j,Ny = sys.lat.Ny, Ny_r = mpi.Int_y[mpi.rank];
    for (i=mpi.start_x; i<mpi.end_x; i++) {
        for (j=0; j<mpi.Int_y[mpi.rank]; j++) {
            receive[ind(j,i,Ny_r)] = send[ind(j+mpi.start_y,i-mpi.start_x,Ny)];
        }
    }
    for (i=0; i<mpi.world_size; i++) {
        if (i != mpi.rank) {
            MPI_Request request; MPI_Status status;
            MPI_Isend(&send[start_index(mpi.Int_y,i)],1,mpi.rows_t[i],i,0,MPI_COMM_WORLD,&request);
            MPI_Recv(&receive[start_index(mpi.Int_x,i)*mpi.Int_y[mpi.rank]],1,mpi.columns_t[i],i,0,MPI_COMM_WORLD,&status);
            MPI_Wait(&request, &status);
        }
    }
}

void scatter_rhs(System sys) {
    MPI mpi = sys.mpi;
    int Nx = sys.lat.Nx, Ny_r = mpi.Int_y[mpi.rank];
    int i, ws = mpi.world_size, rank = mpi.rank;

    if (rank==0) {
        int *counts = malloc(ws * sizeof(int));
        int *displacements = malloc(ws * sizeof(int));
        for (i=0; i<ws; i++) {
            counts[i] = Nx*mpi.Int_y[i];
            displacements[i] = Nx*start_index(mpi.Int_y,i);
        }
        MPI_Scatterv(sys.rhs, counts, displacements, MPI_LONG_DOUBLE, &sys.rhs_r[0], Nx*Ny_r, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
        free(counts); counts = NULL;
        free(displacements); displacements = NULL;
    }
    else {
        MPI_Scatterv(NULL, NULL, NULL, MPI_LONG_DOUBLE, &sys.rhs_r[0], Nx*Ny_r, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void gather_sol(System sys) {
    MPI mpi = sys.mpi;
    int Nx = sys.lat.Nx, Ny_r = mpi.Int_y[mpi.rank];
    int i, ws = mpi.world_size, rank = mpi.rank;
    if (rank==0) {
        int *counts = malloc(ws * sizeof(int));
        int *displacements = malloc(ws * sizeof(int));
        for (i=0; i<ws; i++) {
            counts[i] = Nx*mpi.Int_y[i];
            displacements[i] = Nx*start_index(mpi.Int_y,i);
        }
        MPI_Gatherv(&sys.sol_r[0], Nx*Ny_r, MPI_LONG_DOUBLE, sys.sol, counts, displacements, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
        free(counts); counts = NULL;
        free(displacements); displacements = NULL;
    }
    else {
        MPI_Gatherv(&sys.sol_r[0], Nx*Ny_r, MPI_LONG_DOUBLE, NULL, NULL, NULL, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
    }
}
