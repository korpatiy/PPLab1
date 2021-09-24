#include <iostream>
#include "mpi.h"
#include <cmath>

using namespace std;

const int N = 10000;
const double A = 1;
const double B = 2;
const double H = (B - A) / N;

double func(double x) {
    return 1 / sqrt(2 * x);
}

double calc(int begi, int endi) {
    double result = 0;
    for (int i = begi + 1; i <= endi; i++) {
        double xi = A + H * i;
        result += H * func(xi);
    }
    return result;
}

void get_indexes(int rank, int size, int &begi, int &endi) {
    int p = N / size;
    begi = rank * p;
    if (rank == size - 1)
        endi = N;
    else endi = begi + p;
}

int main(int argc, char **argv) {
    int rank, size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const auto t = MPI_Wtime();
    int begi, endi;
    get_indexes(rank, size, begi, endi);

}
