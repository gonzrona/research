System userInput();
System defineSystem(int argc, char **argv);

double _Complex analytic_solution(double x, double y);
double _Complex f(double x, double y);
double _Complex k(double y);

double _Complex k2(double y);
double _Complex k2_y(double y);
double _Complex k2_yy(double y);
double _Complex f_y(double x, double y);
double _Complex f_xx(double x, double y);
double _Complex f_yy(double x, double y);

void k2_y_4th_order(System sys, double _Complex *k2_e, double _Complex *k2y);
void f_y_4th_order(System sys, double _Complex *f_e, double _Complex *fy);
void f_yy_4th_order(System sys, double _Complex *f_e, double _Complex *fyy_e);
void f_xx_4th_order(System sys, double _Complex *f_e, double _Complex *fxx_e);
void derivative_xx_2nd_order(System sys, double _Complex *f, double _Complex *fxx);
void derivative_yy_2nd_order(System sys, double _Complex *f, double _Complex *fyy);

void coefficients(System sys);
void rhs(System sys);
void LU(System sys);
void solver(System sys);
void residual(System sys);
double normInf(int n, double _Complex *x);

void clearMemory(System sys);

Time tic();
Time toc(Time time);

MPI build_mpi(int argc, char **argv, System sys);
void destroy_mpi(MPI mpi);
int start_index(int *array, int rank);
int end_index(int *array, int rank);
void create_columns(double _Complex *send, double _Complex *receive, System sys);
void create_rows(double _Complex *send, double _Complex *receive, System sys);
void scatter_rhs(System sys);
void gather_sol(System sys);

void printOrder(System sys);
void printLattice(Lattice lat);
void printDoubleArray(int n, double *a, char c[]);
void print_2D_array(double *array, int Nx, int Ny);
void print_2D_array_T(double *array, int Nx, int Ny);
void print_1D_array(double *array, int n, int rank);
void print_1D_array_T(double *array, int Nx, int Ny, int rank);
void print_1D_array_step(double *array, int n, int rank, int step);
