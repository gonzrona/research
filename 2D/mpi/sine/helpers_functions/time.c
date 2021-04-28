#include "../headers/structs.h"

double cpuSecond() {
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ( (double)tp.tv_sec + (double) tp.tv_usec*1.e-6 );
}

Time tic(){
    Time time;
    
    time.time0 = clock();
    time.start_t = MPI_Wtime();
    time.start_t_n = cpuSecond();
    
    return time;
}

Time toc(Time time){
    clock_t time1 = clock();
    double end_t = MPI_Wtime();
    double end_t_n = cpuSecond();

    time.computed_time = (float)(time1 - time.time0)/CLOCKS_PER_SEC;
    time.computed_t = end_t - time.start_t;
    time.computed_t_n = end_t_n - time.start_t_n;
    
    return time;
}


