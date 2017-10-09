#include <iostream>
#include <algorithm>
#include <chrono>

#define N 10000 // N * N = size of matrix

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

    std::cout << "main diagonal - " << "min: "
              << diagMin << "; max: " << diagMax
              << "; avg: " << diagAverage << "\n";
    std::cout << "full matrix - " << "min: " << matrixMin << "; max: " << matrixMax << "\n";
}

int main() {

    int **matrix;
    initMatrix(matrix);

    high_resolution_clock::time_point begin = high_resolution_clock::now();
    calculateMatrixValuesSequentially(matrix);
    high_resolution_clock::time_point end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(end - begin).count();

    std::cout << "Sequential calculation duration: " << duration;

    return 0;
}