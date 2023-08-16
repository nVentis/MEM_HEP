#include <math.h>

double W_jet(double E_j, double E_q, double a[], double b[]){
    double p1 = a[0] + E_q * b[0];
    double p2 = a[1] + E_q * b[1];
    double p3 = a[2] + E_q * b[2];
    double p4 = a[3] + E_q * b[3];
    double p5 = a[4] + E_q * b[4];

    return 1/(
        sqrt(2.*M_PI) * (p2 + p3*p5)
    ) * (exp(
        -pow( (E_j - E_q) - p1 , 2)/(2*pow(p2, 2))
    ) + p3 * exp(
        -pow( (E_j - E_q) - p4 , 2)/(2*pow(p5, 2))
    ));
}

double W_jet(double E_j, double E_q, double JES, double a[], double b[]){
    return (W_jet(E_j/JES, E_q, a, b)/JES);
}