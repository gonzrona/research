#include "../headers/structs.h"

void printLattice(Lattice lat) {
    printf("\n    LATTICE:\n");
    printf("\tNx: %d, Ny: %d\n", lat.Nx, lat.Ny);
    printf("\thx: %f, hy: %f\n", lat.hx, lat.hy);
    printf("\tx0: %f, x1: %f\n", lat.x0, lat.x1);
    printf("\ty0: %f, y1: %f\n\n", lat.y0, lat.y1);
}

void printDoubleArray(int n, double *a, char c[]) {
    int i;
    printf("\n");
    for( i = 0; i < n; i++) {
        printf("\t%s[%d] = %f\n", c, i, a[i]);
    }
    printf("\n");
}

void printOrder(System sys){
    if (sys.order == sixth) {
        printf("\n\t6th Order\n");
    }
    else if (sys.order == fourth) {
        printf("\n\t4th Order\n");
    }
    else {
        printf("\n\t2nd Order\n");
    }
    printf("\tGrid %d x %d\n", sys.lat.Nx, sys.lat.Ny);
}

void print_2D_array(double *array, int Nx, int Ny) {
    int i,j;
    printf("\t[");
    for (j=0; j<Ny; j++) {
        for (i=0; i<Nx; i++) {
            if (array[ind(i,j,Nx)] < 0) {  printf(" %.1f",array[ind(i,j,Nx)]); }
            else { printf("  %.1f",array[ind(i,j,Nx)]); }
        }
        if (j != Ny-1) { printf("\n\t "); }
    }
    printf(" ]\n");
}
void print_1D_array(double *array, int n, int rank) {
    int i;
    printf("\trank %d: [",rank);
    for (i=0; i<n; i++) {
        if (array[i] < 0) {  printf(" %.1f",array[i]); }
        else { printf("  %.1f",array[i]); }
    }
    printf(" ]\n");
}
void print_1D_array_step(double *array, int n, int rank, int step) {
    int i;
    printf("\tstep %d, rank %d: [",step,rank);
    for (i=0; i<n; i++) {
        if (array[i] < 0) {  printf(" %.1f",array[i]); }
        else { printf("  %.1f",array[i]); }
    }
    printf(" ]\n");
}
void print_2D_array_T(double *array, int Nx, int Ny) {
    int i,j;
    printf("\t[");
    for (j=0; j<Ny; j++) {
        for (i=0; i<Nx; i++) {
            if (array[ind(j,i,Ny)] < 0) {  printf(" %.1f",array[ind(j,i,Ny)]); }
            else { printf("  %.1f",array[ind(j,i,Ny)]); }
        }
        if (j != Ny-1) { printf("\n\t "); }
    }
    printf(" ]'\n");
}
void print_1D_array_T(double *array, int Nx, int Ny, int rank) {
    int i,j;
    printf("\trank %d: [",rank);
    for (j=0; j<Ny; j++) {
        for (i=0; i<Nx; i++) {
            if (array[ind(j,i,Ny)] < 0) {  printf(" %.1f",array[ind(j,i,Ny)]); }
            else { printf("  %.1f",array[ind(j,i,Ny)]); }
        }
    }
    printf(" ]'\n");
}
