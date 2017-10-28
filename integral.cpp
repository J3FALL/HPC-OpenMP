#include "integral.h"
#include <math.h>
#include <iostream>
#include <omp.h>


double fun(double x) {
    return pow(sin(1.0 / x), 2.0) / pow(x, 2.0);
}

bool isEqualOrder(double a, double b) {
    return (a / b < 10.0 && a / b > 1.0) ? true : false;
}

double calculateSerialIntegral(double start, double end, double precision) {

    //initial approximation for 2 steps
    long step = 2;
    double prev_approx = 0.0;
    double current_approx = ((end - start) / step) * (fun(start) / step + fun(end) / step);

    while (!isEqualOrder(fabs(prev_approx - current_approx), precision)) {
        std::cout << fabs(prev_approx - current_approx) << " " << precision << "\n";
        std::cout << "Current steps: " << step << "\n";
        prev_approx = current_approx;
        step *= 2;

        double middle_sum = 0.0;
        for (long i = 1; i <= step - 1; i++) {
            double x = start + (end - start) * i / step;
            middle_sum += fun(x);
        }

        current_approx = ((end - start) / step) * (fun(start) / step + middle_sum + fun(end) / step);
    }

    std::cout << "Number of steps: " << step << "\n";
    return current_approx;
}

double calculateIntegralWithAtomic(double start, double end, double precision) {

    int numberOfThreads = 8;

    omp_set_num_threads(numberOfThreads);

    long step = 2;
    double prev_approx = 0.0;
    double current_approx = ((end - start) / step) * (fun(start) / step + fun(end) / step);

    while (!isEqualOrder(fabs(prev_approx - current_approx), precision)) {
        std::cout << fabs(prev_approx - current_approx) << " " << precision << "\n";
        std::cout << "Current steps: " << step << "\n";
        prev_approx = current_approx;
        step *= 2;

        double middle_sum = 0.0;

#pragma omp parallel for schedule(static)
        for (long i = 1; i <= step - 1; i++) {
            double x = start + (end - start) * i / step;

#pragma omp atomic
            middle_sum += fun(x);
        }

        current_approx = ((end - start) / step) * (fun(start) / step + middle_sum + fun(end) / step);
    }

    std::cout << "Number of steps: " << step << "\n";
    return current_approx;
}

double calculateIntegralWithCritical(double start, double end, double precision) {

    int numberOfThreads = 8;

    omp_set_num_threads(numberOfThreads);

    long step = 2;
    double prev_approx = 0.0;
    double current_approx = ((end - start) / step) * (fun(start) / step + fun(end) / step);

    while (!isEqualOrder(fabs(prev_approx - current_approx), precision)) {
        std::cout << fabs(prev_approx - current_approx) << " " << precision << "\n";
        std::cout << "Current steps: " << step << "\n";
        prev_approx = current_approx;
        step *= 2;

        double middle_sum = 0.0;

#pragma omp parallel for schedule(static)
        for (long i = 1; i <= step - 1; i++) {
            double x = start + (end - start) * i / step;

#pragma omp critical
            middle_sum += fun(x);
        }

        current_approx = ((end - start) / step) * (fun(start) / step + middle_sum + fun(end) / step);
    }

    std::cout << "Number of steps: " << step << "\n";
    return current_approx;
}

double calculateIntegralWithLocks(double start, double end, double precision) {

    int numberOfThreads = 8;

    omp_set_num_threads(numberOfThreads);

    long step = 2;
    double prev_approx = 0.0;
    double current_approx = ((end - start) / step) * (fun(start) / step + fun(end) / step);

    omp_lock_t sum_lock;
    omp_init_lock(&sum_lock);

    while (!isEqualOrder(fabs(prev_approx - current_approx), precision)) {
        std::cout << fabs(prev_approx - current_approx) << " " << precision << "\n";
        std::cout << "Current steps: " << step << "\n";
        prev_approx = current_approx;
        step *= 2;

        double middle_sum = 0.0;

#pragma omp parallel for schedule(static)
        for (long i = 1; i <= step - 1; i++) {
            double x = start + (end - start) * i / step;

            omp_set_lock(&sum_lock);

            middle_sum += fun(x);

            omp_unset_lock(&sum_lock);
        }

        current_approx = ((end - start) / step) * (fun(start) / step + middle_sum + fun(end) / step);
    }

    omp_destroy_lock(&sum_lock);

    std::cout << "Number of steps: " << step << "\n";
    return current_approx;
}

double calculateIntegralWithReduction(double start, double end, double precision) {

    int numberOfThreads = 8;

    omp_set_num_threads(numberOfThreads);

    long step = 2;
    double prev_approx = 0.0;
    double current_approx = ((end - start) / step) * (fun(start) / step + fun(end) / step);

    while (!isEqualOrder(fabs(prev_approx - current_approx), precision)) {
        std::cout << fabs(prev_approx - current_approx) << " " << precision << "\n";
        std::cout << "Current steps: " << step << "\n";
        prev_approx = current_approx;
        step *= 2;

        double middle_sum = 0.0;

#pragma omp parallel for schedule(static) reduction(+:middle_sum)
        for (long i = 1; i <= step - 1; i++) {
            double x = start + (end - start) * i / step;
            middle_sum += fun(x);
        }

        current_approx = ((end - start) / step) * (fun(start) / step + middle_sum + fun(end) / step);
    }

    std::cout << "Number of steps: " << step << "\n";
    return current_approx;
}
