#include "../headers/structs.h"
#include "../headers/prototypes.h"

System defineSystem(int argc, char **argv) {
    
    System sys = userInput();

    if (argc == 3) {
        if (atoi(argv[1]) == 6) { sys.order = sixth; }
        else if (atoi(argv[1]) == 4) { sys.order = fourth; }
        else { sys.order = second; }
        sys.lat.Nx = atoi(argv[2]);
        sys.lat.Ny = sys.lat.Nx;
    }
    else if (argc == 2) {
        if (atoi(argv[1]) == 6) { sys.order = sixth; }
        else if (atoi(argv[1]) == 4) { sys.order = fourth; }
        else { sys.order = second; }
        sys.lat.Nx = 1100;
        sys.lat.Ny = 1000;
    }
    else {
        sys.order = second;
        sys.lat.Nx = 1100;
        sys.lat.Ny = 1000;
    }
    
    if (sys.order == sixth && sys.lat.Nx != sys.lat.Ny) { sys.lat.Ny = sys.lat.Nx; }
    
    sys.mpi = build_mpi(argc, argv, sys);

    if (sys.order == sixth && sys.lat.Nx != sys.lat.Ny) {
        if (sys.mpi.rank==0){printf("WARNING: Sixth order requires a uniform grid size, Ny was set to Nx\n");}
    }

    sys.lat.hx = (sys.lat.x1-sys.lat.x0)/(sys.lat.Nx+1);
    sys.lat.hy = (sys.lat.y1-sys.lat.y0)/(sys.lat.Ny+1);
    
    sys.lat.Nxy = sys.lat.Nx * sys.lat.Ny;

    sys.sol = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex));
    sys.sol_analytic = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex));
    sys.rhs = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex));
    sys.L = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex));
    sys.U = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex));
    sys.Up = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex));
    sys.res = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex));
    sys.error = (double _Complex *) malloc(sys.lat.Nxy * sizeof(double _Complex)); // can reuse res rather than allocating more memory
    sys.k = (double _Complex *) malloc(sys.lat.Ny * sizeof(double _Complex));
    
    sys.bp = (double _Complex *) malloc(sys.lat.Ny * sizeof(double _Complex));
    sys.ap = (double _Complex *) malloc(sys.lat.Ny * sizeof(double _Complex));
    sys.b  = (double _Complex *) malloc(sys.lat.Ny * sizeof(double _Complex));
    sys.a  = (double _Complex *) malloc(sys.lat.Ny * sizeof(double _Complex));
    sys.bm = (double _Complex *) malloc(sys.lat.Ny * sizeof(double _Complex));
    sys.am = (double _Complex *) malloc(sys.lat.Ny * sizeof(double _Complex));
    
    sys.rhs_r = (double _Complex *) malloc(sys.lat.Nx*sys.mpi.Int_y[sys.mpi.rank] * sizeof(double _Complex));
    sys.sol_r = (double _Complex *) malloc(sys.lat.Nx*sys.mpi.Int_y[sys.mpi.rank] * sizeof(double _Complex));

    
    
    return sys;
}

void clearMemory(System sys) {
    free(sys.sol); sys.sol = NULL;
    free(sys.sol_analytic); sys.sol_analytic = NULL;
    free(sys.rhs); sys.rhs = NULL;
    free(sys.L); sys.L = NULL;
    free(sys.U); sys.U = NULL;
    free(sys.Up); sys.Up = NULL;
    free(sys.res); sys.res = NULL;
    free(sys.error); sys.error = NULL; // can reuse res rather than allocating more memory
    free(sys.k); sys.k = NULL;
    
    free(sys.bp); sys.bp = NULL;
    free(sys.ap); sys.ap = NULL;
    free(sys.b);  sys.b  = NULL;
    free(sys.a);  sys.a  = NULL;
    free(sys.bm); sys.bm = NULL;
    free(sys.am); sys.am = NULL;

    free(sys.rhs_r); sys.rhs_r = NULL;
    free(sys.sol_r); sys.sol_r = NULL;

}
