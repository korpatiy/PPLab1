#include <iostream>
#include <cmath>
#include "mpi.h"

using namespace std;

const int N = 1000;
const double A = 1;
const double B = 2;
const double H = (B - A) / N;

double func(double x) {
    return 1 / sqrt(2 * x);
}

double calc(int begi, int endi) {
    double result = 0;
    for (int i = begi; i < endi; i++) {
        double xi = A + H * i;
        result += H * func(xi);
    }
    return result;
}

void get_indexes(int rank, int size, int& begi, int& endi) {
    int p = N / size;
    begi = rank * p;
    if (rank == size - 1)
        endi = N;
    else endi = begi + p;
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const auto t = MPI_Wtime();
    int begi, endi;
    get_indexes(rank, size, begi, endi);
    const auto result = calc(begi, endi);

    if (rank == 0) {
        double sum = result;
        for (int i = 1; i < size; i++) {
            double part_sum = 0;
            MPI_Recv(&part_sum, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
            sum += part_sum;
        }
        cout << "Sum = " << sum << endl;
    }
    else {
        MPI_Send(&result, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    cout << "Rank = " << rank << ": time = " << (MPI_Wtime() - t) << " sec" << endl;
    MPI_Finalize();
    clock_t t1 = clock();
    const auto result1 = calc(0, N);
    cout << "Sum = " << result1 << endl;
    double time = (static_cast<double>(clock()) - t1) / CLOCKS_PER_SEC;
    cout << "Time = " << time << " sec" << endl;
}
