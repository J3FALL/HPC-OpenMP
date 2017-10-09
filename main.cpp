#include <iostream>
#include <time.h>

#define N 10000 // N * N = size of matrix

void initMatrix(int **&matrix) {

    srand(time(NULL));

    matrix = new int *[N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new int[N];
        for (int j = 0; j < N; j++) {
            matrix[i][j] = rand() % 10000;
        }
    }
}

int main() {

    int **matrix;
    initMatrix(matrix);

}