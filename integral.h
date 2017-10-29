

#ifndef HPC_OPENMP_INTEGRAL_H
#define HPC_OPENMP_INTEGRAL_H

struct IntegralResults {
    double result;
    double points;
};

IntegralResults calculateSerialIntegral(double start, double end, double precision);

IntegralResults calculateIntegralWithAtomic(double start, double end, double precision);

IntegralResults calculateIntegralWithCritical(double start, double end, double precision);

IntegralResults calculateIntegralWithLocks(double start, double end, double precision);

IntegralResults calculateIntegralWithReduction(double start, double end, double precision);

#endif //HPC_OPENMP_INTEGRAL_H
