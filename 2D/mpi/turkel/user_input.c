#include "headers/structs.h"
#include "headers/prototypes.h"

Parameters p;

double _Complex k(double y) {
    double a = p.a, b = p.b, c = p.c;
    return a - b * sin(c * y);
}

double _Complex k2(double y) {
    return k(y)*k(y);
}

double _Complex k2_y(double y) {
    double b = p.b, c = p.c;
    return -2.0*b*c*k(y)*cos(c*y);
}

double _Complex k2_yy(double y) {
    double b = p.b, c = p.c;
    return 2.0*b*c*c*(k(y)*sin(c*y) + b*cos(c*y)*cos(c*y));
}

double _Complex analytic_solution(double x, double y) {
    double c = p.c, beta = p.beta;
    return sin(beta * x) * exp(-k(y)/c);
}

double _Complex f(double x, double y) {
    double a = p.a, b = p.b, c = p.c, beta = p.beta;
    return -b * (2.0*a + c) * sin(c*y) * exp(-k(y)/c) * sin(beta * x);
}

double _Complex f_xx(double x, double y) {
    double a = p.a, b = p.b, c = p.c, beta = p.beta;
    return b * (2.0*a + c) * beta*beta * sin(c*y) * exp(-k(y)/c) * sin(beta*x);
}

double _Complex f_y(double x, double y) {
    double a = p.a, b = p.b, c = p.c, beta = p.beta;
    return -b * (2.0*a + c) * sin(beta*x) * ( b * sin(c*y) + c ) * exp(-k(y)/c) * cos(c*y);
}

double _Complex f_yy(double x, double y) {
    double a = p.a, b = p.b, c = p.c, beta = p.beta;
    double sin2 = sin(c*y) * sin(c*y), cos2 = cos(c*y) * cos(c*y);
    return -b*(2.0*a+c)*sin(beta*x) * exp(-k(y)/c) * ( -b*c*sin2 - c*c*sin(c*y) + b*b*sin(c*y)*cos2 + 2.0*b*c*cos2 );
}


System userInput() {
    System sys;

//    DEFINE THE DOMAIN
    sys.lat.x0 = 0.0; sys.lat.x1 = M_PI;
    sys.lat.y0 = 0.0; sys.lat.y1 = M_PI;

    
    //    DEFINE PARAMETERS
    p.a = 10.0;
    p.b = 9.0;
    p.c = 10.0;
    p.beta = sqrt(p.a * p.a + p.b * p.b);

    return sys;
}
