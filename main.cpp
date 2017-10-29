#include <iostream>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include "integral.h"

#define N 10000 // N * N = size of matrix
#define VECTOR_SIZE 100000000

int numberOfThreads = 16;

using namespace std::chrono;

void initMatrix(int **&matrix) {

    srand(time(NULL));

    matrix = new int *[N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new int[N];
        for (int j = 0; j < N; j++) {
            matrix[i][j] = -250 + rand() % 1000;
        }
    }
}

void calculateMatrixValuesSequentially(int **matrix) {
    std::cout << "Starting sequential calculations \n";

    int rowsMin[N], rowsMax[N];
    double rowsAverage[N];
    int diagMin = matrix[0][0], diagMax = matrix[0][0];
    int diagSum = 0;
    double diagAverage;
    int matrixMax, matrixMin;

    for (int i = 0; i < N; i++) {

        int rowSum = 0;
        int rowMin = matrix[i][0];
        int rowMax = matrix[i][0];

        for (int j = 0; j < N; j++) {
            rowSum += matrix[i][j];

            rowMin = std::min(rowMin, matrix[i][j]);
            rowMax = std::max(rowMax, matrix[i][j]);

            if (i == j) {
                //main diagonal
                diagMin = std::min(diagMin, matrix[i][j]);
                diagMax = std::max(diagMax, matrix[i][j]);
                diagSum += matrix[i][j];
            }
        }

        double avg = (double) rowSum / N;
        rowsAverage[i] = avg;
        rowsMin[i] = rowMin;
        rowsMax[i] = rowMax;

        matrixMin = std::min(matrixMin, rowMin);
        matrixMax = std::max(matrixMax, rowMax);

    }

    diagAverage = (double) diagSum / N;

    std::cout << "main diagonal sum: " << diagSum << "\n";
    std::cout << "main diagonal - " << "min: "
              << diagMin << "; max: " << diagMax
              << "; avg: " << diagAverage << "\n";
    std::cout << "full matrix - " << "min: " << matrixMin << "; max: " << matrixMax << "\n";
}

void calculateMatrixValuesParallel(int **matrix) {

    omp_set_num_threads(numberOfThreads);

    std::cout << "Starting parallel calculations \n";

    int rowsMin[N], rowsMax[N];
    double rowsAverage[N];
    int diagMin = matrix[0][0], diagMax = matrix[0][0];
    long diagSum = 0;
    double diagAverage;
    int matrixMax, matrixMin;

#pragma omp parallel for schedule(static)  reduction(+:diagSum)
    for (int i = 0; i < N; i++) {

        //std::cout << "iteration: " << i << " from thread: " << omp_get_thread_num() << "\n";
        int rowSum = 0;
        int rowMin, rowMax;

        rowMin = matrix[i][0];
        rowMax = matrix[i][0];

        for (int j = 0; j < N; j++) {
            rowSum += matrix[i][j];
            rowMin = std::min(rowMin, matrix[i][j]);
            rowMax = std::max(rowMax, matrix[i][j]);

            if (i == j) {
                //main diagonal
                diagMin = std::min(diagMin, matrix[i][j]);
                diagMax = std::max(diagMax, matrix[i][j]);

                diagSum += matrix[i][j];

            }
        }

        double avg = (double) rowSum / N;

        rowsAverage[i] = avg;
        rowsMin[i] = rowMin;
        rowsMax[i] = rowMax;

        matrixMin = std::min(matrixMin, rowMin);
        matrixMax = std::max(matrixMax, rowMax);

    }

    diagAverage = (double) diagSum / N;

    std::cout << "main diagonal sum: " << diagSum << "\n";
    std::cout << "main diagonal - " << "min: "
              << diagMin << "; max: " << diagMax
              << "; avg: " << diagAverage << "\n";
    std::cout << "full matrix - " << "min: " << matrixMin << "; max: " << matrixMax << "\n";
}

void runMatrixExperiment() {
    std::cout << "Run matrix experiment\n";

    int **matrix;
    initMatrix(matrix);

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    calculateMatrixValuesSequentially(matrix);
    high_resolution_clock::time_point end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Sequential calculation duration: " << duration << "\n";

    begin = high_resolution_clock::now();
    calculateMatrixValuesParallel(matrix);
    end = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Parallel calculation duration: " << duration << "\n";
}

long calculateScalarProductSequentially(int *vectorA, int *vectorB) {
    std::cout << "Starting sequential calculations \n";

    long product = 0;

    for (int i = 0; i < VECTOR_SIZE; i++) {
        product += vectorA[i] * vectorB[i];
    }

    return product;
}

long calculateScalarProductParallel(int *vectorA, int *vectorB) {
    std::cout << "Starting parallel calculations \n";

    omp_set_num_threads(numberOfThreads);

    long product = 0;

#pragma omp parallel for schedule(static) reduction(+:product)
    for (int i = 0; i < VECTOR_SIZE; i++) {
        product += vectorA[i] * vectorB[i];
    }

    return product;
}

void initVector(int *&vector) {

    srand(time(NULL));

    vector = new int[VECTOR_SIZE];
    for (int i = 0; i < VECTOR_SIZE; i++) {
        vector[i] = -250 + rand() % 1000;
    }
}

void runVectorExperiment() {
    std::cout << "Run vector experiment\n";

    int *vectorA, *vectorB;

    initVector(vectorA);
    initVector(vectorB);

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    long product = calculateScalarProductSequentially(vectorA, vectorB);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Product: " << product << "\n";
    std::cout << "Sequential calculation duration: " << duration << "\n";

    begin = high_resolution_clock::now();
    product = calculateScalarProductParallel(vectorA, vectorB);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Product: " << product << "\n";
    std::cout << "Parallel calculation duration: " << duration << "\n";
}

void runSerialIntegralExperiment() {
    std::cout << "Run serial integral calculation experiment\n";

    double eps = 1.0e-06;
    double a = 10.0, b = 100.0;

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    double result = calculateSerialIntegral(a, b, eps);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Result: " << result << "\n";
    std::cout << "Serial integral calculation duration: " << duration << "\n";

}

void runIntegralWithAtomicExperiment() {
    std::cout << "Run integral calculation with atomic experiment\n";

    double eps = 1.0e-06;
    double a = 0.00001, b = 0.0001;

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    double result = calculateIntegralWithAtomic(a, b, eps);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Result: " << result << "\n";
    std::cout << "Integral calculation with atomic duration: " << duration << "\n";

}

void runIntegralWithCriticalExperiment() {
    std::cout << "Run integral calculation with critical experiment\n";

    double eps = 1.0e-06;
    double a = 0.00001, b = 0.0001;

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    double result = calculateIntegralWithCritical(a, b, eps);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Result: " << result << "\n";
    std::cout << "Integral calculation with critical duration: " << duration << "\n";

}

void runIntegralWithLocksExperiment() {
    std::cout << "Run integral calculation with locks experiment\n";

    double eps = 1.0e-06;
    double a = 0.00001, b = 0.0001;

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    double result = calculateIntegralWithCritical(a, b, eps);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Result: " << result << "\n";
    std::cout << "Integral calculation with locks duration: " << duration << "\n";

}

void runIntegralWithReductionExperiment() {
    std::cout << "Run integral calculation with reduction experiment\n";

    double eps = 1.0e-06;
    double a = 0.00001, b = 0.0001;

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    double result = calculateIntegralWithCritical(a, b, eps);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Result: " << result << "\n";
    std::cout << "Integral calculation with reduction duration: " << duration << "\n";

}

int main() {
    /*
    runMatrixExperiment();
    runVectorExperiment();
    */

    runSerialIntegralExperiment();
    runIntegralWithAtomicExperiment();
    runIntegralWithCriticalExperiment();
    runIntegralWithLocksExperiment();
    runIntegralWithReductionExperiment();

    return 0;
}