#include <iostream>
#include "mpi.h"
#include <stdio.h>

#define N 256

void print_results(char* prompt, int result[N][N]);

int main(int argc, char* argv[])
{
    int rank; // Ранг процесса
    int size; // Номер процесса
    int i, j, k, rows = 0;
    int buff = 0;

    int mA[N][N];
    int mB[N][N];
    int mResult[N][N];

    double startTime, endTime;

    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv); //initialize MPI operations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) 
    {
        // Инициализация матриц А и B
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                mA[i][j] = rand() % 10;
        
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                mB[i][j] = rand() % 10;

        startTime = MPI_Wtime();
        for (i = 0; i < size; i++)
        {
            if (i < (N % (size - 1)))
                rows = N / (size - 1) + 1;
            else 
                rows = N / (size - 1);

            MPI_Send(&buff, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&mA[buff][0], N * N, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&mB, N * N, MPI_INT, i, 1, MPI_COMM_WORLD);
            buff += rows;
        }

        for (i = 1; i < size; i++) {
            MPI_Recv(&buff, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&mResult[buff][0], rows * N, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
        }
        endTime = MPI_Wtime();
        printf("\nRunning Time = %f\n\n", endTime - startTime);
    }

    if (rank > 0)
    {
        MPI_Recv(&buff, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&mA[buff][0], N * N, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&mB, N * N, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        for (i = buff; i < rows; i++)
            for (k = 0; k < N; k++)
                for (j = 0; j < N; j++)
                    mResult[i][j] += (mA[i][k] * mB[k][j]);

        MPI_Send(&buff, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&(mResult[buff][0]), rows * N, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }


    MPI_Finalize();
    return 0;

}