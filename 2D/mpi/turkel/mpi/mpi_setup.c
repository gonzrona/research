#include "../headers/structs.h"

int start_index(int *array, int rank) {
    int start = 0;
    int i;
    for (i=0; i<rank; i++) {
        start = start + array[i];
    }
    return start;
}
int end_index(int *array, int rank) {
    int end = 0;
    int i;
    for (i=0; i<=rank; i++) {
        end = end + array[i];
    }
    return end;
}

MPI build_mpi(int argc, char **argv, System sys) {
    MPI mpi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.world_size);
    
    int i, Nx = sys.lat.Nx, Ny = sys.lat.Ny, world_size = mpi.world_size;

    mpi.Int_x = (int*) malloc(world_size * sizeof(int));
    mpi.Int_y = (int*) malloc(world_size * sizeof(int));
    mpi.columns = (MPI_Datatype*) malloc(world_size * sizeof(MPI_Datatype));
    mpi.rows = (MPI_Datatype*) malloc(world_size * sizeof(MPI_Datatype));
    mpi.columns_t = (MPI_Datatype*) malloc(world_size * sizeof(MPI_Datatype));
    mpi.rows_t = (MPI_Datatype*) malloc(world_size * sizeof(MPI_Datatype));

    for (i=0; i<Nx%world_size; i++) {
        mpi.Int_x[i]= Nx/world_size+1;
    }
    for (i=Nx%world_size; i<world_size; i++) {
        mpi.Int_x[i]= Nx/world_size;
    }
    for (i=0; i<Ny%world_size; i++) {
        mpi.Int_y[i]= Ny/world_size+1;
    }
    for (i=Ny%world_size; i<world_size; i++) {
        mpi.Int_y[i]= Ny/world_size;
    }

    mpi.start_y = start_index(mpi.Int_y,mpi.rank);
    mpi.end_y = end_index(mpi.Int_y,mpi.rank);
    mpi.start_x = start_index(mpi.Int_x,mpi.rank);
    mpi.end_x = end_index(mpi.Int_x,mpi.rank);
    
    for (i=0; i<world_size; i++) {
        if (i != mpi.rank) {
            MPI_Datatype column;
            MPI_Type_vector(mpi.Int_y[mpi.rank], mpi.Int_x[i], Nx, MPI_LONG_DOUBLE, &column);
            MPI_Type_commit(&column);
            mpi.columns[i] = column;

            MPI_Datatype row;
            MPI_Type_contiguous(mpi.Int_y[i]*mpi.Int_x[mpi.rank], MPI_LONG_DOUBLE, &row);
            MPI_Type_commit(&row);
            mpi.rows[i] = row;

            MPI_Datatype column_t;
            MPI_Type_contiguous(mpi.Int_x[i]*mpi.Int_y[mpi.rank], MPI_LONG_DOUBLE, &column_t);
            MPI_Type_commit(&column_t);
            mpi.columns_t[i] = column_t;

            MPI_Datatype row_t;
            MPI_Type_vector(mpi.Int_x[mpi.rank], mpi.Int_y[i] , Ny, MPI_LONG_DOUBLE, &row_t);
            MPI_Type_commit(&row_t);
            mpi.rows_t[i] = row_t;
          }
    }
    return mpi;
}

void destroy_mpi(MPI mpi) {
    int i;
    for (i=0; i<mpi.world_size; i++) {
        if (i != mpi.rank) {
            MPI_Type_free (&mpi.columns[i]);
            MPI_Type_free (&mpi.columns_t[i]);
            MPI_Type_free (&mpi.rows[i]);
            MPI_Type_free (&mpi.rows_t[i]);
        }
    }
    free(mpi.Int_x); mpi.Int_x = NULL;
    free(mpi.Int_y); mpi.Int_y = NULL;
    free(mpi.columns); mpi.columns = NULL;
    free(mpi.rows); mpi.rows = NULL;
    free(mpi.columns_t); mpi.columns_t = NULL;
    free(mpi.rows_t); mpi.rows_t = NULL;
    MPI_Finalize();
}

