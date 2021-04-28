#include "../headers/structs.h"
#include "../headers/prototypes.h"

void DST(DSTN dst, double _Complex *b, double _Complex *bhat, fftw_plan plan, double *in, fftw_complex *out);

void solver(System sys) {
    
    DSTN dst;
    MPI mpi = sys.mpi;
    int i,j,my,mx;
    int Nx = sys.lat.Nx, Ny = sys.lat.Ny, Nx_r = mpi.Int_x[mpi.rank], Ny_r = mpi.Int_y[mpi.rank];
    int start_x = mpi.start_x, end_x = mpi.end_x;
    int N = 2*Nx + 2, NC = (N/2) + 1;
    dst.Nx = Nx; dst.N = N; dst.coef = sqrt(2.0/(Nx+1));

    double _Complex *rhat_h = (double _Complex *) malloc(Nx*Ny_r * sizeof(double _Complex));
    double _Complex *rhat_v = (double _Complex *) malloc(Nx_r*Ny * sizeof(double _Complex));
    double _Complex *xhat_v = (double _Complex *) malloc(Nx_r*Ny * sizeof(double _Complex));
    double _Complex *b      = (double _Complex *) malloc(Nx * sizeof(double _Complex));
    double _Complex *bhat   = (double _Complex *) malloc(Nx * sizeof(double _Complex));
    double _Complex *y    = (double _Complex *) malloc(Ny  * sizeof(double _Complex));
    double *in    = (double *) fftw_malloc(sizeof(double) * N);

    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NC);
    fftw_plan plan  = fftw_plan_dft_r2c_1d ( N, in, out, FFTW_ESTIMATE );
        
    for(j = 0; j < Ny_r; j++) {
        my = j*Nx ;
        for(i = 0; i < Nx; i++){
            b[i] = sys.rhs_r[i + my] ;
        }
        DST(dst, b, bhat, plan, in, out);
        for(i = 0; i < Nx; i++){
            rhat_h[i + my] = bhat[i] ;
        }
    }
    
    create_columns(rhat_h,rhat_v,sys);
    
    for(i = start_x; i < end_x; i++){
        y[0] = rhat_v[i-start_x];
        mx = i*Ny ;
        for(j = 1; j < Ny; j++) {
            y[j] = rhat_v[i-start_x + j*Nx_r] - sys.L[j + mx]*y[j - 1];
        }
        xhat_v[Ny - 1 + (i-start_x)*Ny] = y[Ny - 1]/sys.U[Ny - 1 + mx] ;
        for(j = Ny-2; j >= 0; j--) {
            xhat_v[j + (i-start_x)*Ny] =  ( y[j] - sys.Up[j + mx] * xhat_v[j + 1 + (i-start_x)*Ny] )/sys.U[j + mx] ;
        }
    }
    
    create_rows(xhat_v,rhat_h,sys);
    
    for(j = 0; j < Ny_r; j++) {
        my = j*Nx;
        for(i = 0; i < Nx; i++) {
            b[i] = rhat_h[j + i*Ny_r];
        }
        DST(dst, b, bhat, plan, in, out);
        for(i = 0; i < Nx; i++){
            sys.sol_r[i + my] = bhat[i];
        }
    }
        
    fftw_destroy_plan(plan);
    fftw_free(in); in = NULL;
    fftw_free(out); out = NULL;
    free(b); b = NULL;
    free(bhat); bhat = NULL;
    free(y); y = NULL;
        
    free(rhat_v); rhat_v = NULL;
    free(xhat_v); xhat_v = NULL;
    free(rhat_h); rhat_h = NULL;
}
