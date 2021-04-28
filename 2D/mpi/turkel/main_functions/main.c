#include "../headers/structs.h"
#include "../headers/prototypes.h"

int main (int argc, char **argv) {
    
    System sys = defineSystem(argc, argv);
    MPI mpi = sys.mpi;
    int rank = mpi.rank;
    
    if (rank == 0) { printOrder(sys); }
        
    Time time = tic();

    coefficients(sys);
        
    rhs(sys);

    LU(sys);
    
    time = toc(time);
    
    if (rank == 0) {
        printf("\n \tSet-up time = %f sec \n", time.computed_time );
        printf("   \tSet-up wall-time = %f sec \n", time.computed_t);
        printf("   \tSet-up wall-time (sys/time) = %f sec \n\n", time.computed_t_n);
    }
        
    scatter_rhs(sys);

    time = tic();
        
    solver(sys);

    time = toc(time);
    
    gather_sol(sys);
    
    if (rank == 0) {
        printf("\n \tSolver time = %f sec \n", time.computed_time );
        printf("   \tSolver wall-time = %f sec \n", time.computed_t);
        printf("   \tSolver wall-time (sys/time) = %f sec \n\n", time.computed_t_n);

        residual(sys);

        printf("\t||res||_inf =  %10.7e \n\n",normInf(sys.lat.Nxy, sys.res));
        printf("\t||error||_inf =  %10.7e \n\n",normInf(sys.lat.Nxy, sys.error));
    }

    clearMemory(sys);
    destroy_mpi(mpi);

    return 0;
}
