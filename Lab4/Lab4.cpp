﻿#include <iostream>
#include "mpi.h"
#include <stdio.h>

#define ROW 10
#define COL 10

double mA[ROW][COL];
double mB[ROW][COL];
double mR[ROW][COL];

double startTime;
double endTime;

void MakeMatrix(double matrix[ROW][COL], int n)
{
    for (int i = 0; i < ROW; i++) {
        for (int j = 0; j < COL; j++) {
            matrix[i][j] = i + j * n;
        }
    }
}
void printArray()
{
    int i, j;
    for (i = 0; i < ROW; i++) {
        printf("\n");
        for (j = 0; j < COL; j++)
            printf("%8.2f  ", mA[i][j]);
    }
    printf("\n---------------\n");
    for (i = 0; i < ROW; i++) {
        printf("\n");
        for (j = 0; j < COL; j++)
            printf("%8.2f  ", mB[i][j]);
    }
    printf("\n---------------\n");
    for (i = 0; i < ROW; i++) {
        printf("\n");
        for (j = 0; j < COL; j++)
            printf("%8.2f  ", mR[i][j]);
    }
    printf("\n\n");
}

int main(int argc, char* argv[])
{
    int rank, size;
    int numCalculate; // Кол-во вычислений на один процесс
    int lowBound, upBound;
    int i;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0)
    {
        MakeMatrix(mA, 1);
        MakeMatrix(mB, 2);
        startTime = MPI_Wtime();
        numCalculate = ROW / (size - 1);
        printf("OK");
        lowBound = 0;
        for (i = 1; i < size; i++)
        {
            lowBound = (i - 1) * numCalculate;
            if (((i + 1) == size) && ((ROW % (size - 1)) != 0)) upBound = ROW;
            else upBound = lowBound + numCalculate;
            MPI_Send(&lowBound, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&upBound, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&mA[lowBound][0], (upBound-lowBound) * COL, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send(&mB, ROW * COL, MPI_INT, i, 1, MPI_COMM_WORLD);
        }
    }
    if (rank > 0)
    {
        MPI_Recv(&lowBound, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&upBound, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&mA[lowBound][0], (upBound - lowBound) * COL, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&mB, ROW * COL, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
        printf("rank: %d \tsize: %d - [%d;%d]", rank, size, lowBound, upBound);
        for (i = lowBound; i < upBound; i++)
            for (int k = 0; k < ROW; k++)
                for (int j = 0; j < COL; j++) 
                    mR[i][j] += (mA[i][k] * mB[k][j]);
        MPI_Send(&lowBound, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&upBound, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&mR[lowBound][0], (upBound - lowBound) * COL, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) 
    {
        for (i = 1; i < size; i++)
        {
            MPI_Recv(&lowBound, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&upBound, 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&mR[lowBound][0], (upBound - lowBound) * COL, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
        }
        endTime = MPI_Wtime() - startTime;
        printf("\nTime Parallel = %f\n\n", endTime);
        printArray();
    }
    MPI_Finalize();
    return 0;
}

