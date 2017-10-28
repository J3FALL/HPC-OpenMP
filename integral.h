

#ifndef HPC_OPENMP_INTEGRAL_H
#define HPC_OPENMP_INTEGRAL_H

double calculateSerialIntegral(double start, double end, double precision);

double calculateIntegralWithAtomic(double start, double end, double precision);

double calculateIntegralWithCritical(double start, double end, double precision);

double calculateIntegralWithLocks(double start, double end, double precision);

double calculateIntegralWithReduction(double start, double end, double precision);

#endif //HPC_OPENMP_INTEGRAL_H
