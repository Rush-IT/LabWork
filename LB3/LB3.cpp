#include <iostream>
#include <omp.h>
int main()
{
    // Последовательная реализация метода Гауса
    int n = 1000;
    double t = omp_get_wtime();
    double* a = (double*)malloc(sizeof(*a) * n * n);
    double* b = (double*)malloc(sizeof(*b) * n);
    double* x = (double*)malloc(sizeof(*x) * n);
    
    for (int i = 0; i < n; i++) {
        srand(i * (n + 1));
        for (int j = 0; j < n; j++)
            a[i * n + j] = rand() % 100 + 1;
        b[i] = rand() % 100 + 1;
    }
    #if 0
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%12.4f ", a[i * n + j]);
        printf(" | %12.4f\n ", b[i]);
    }
    #endif

    // Прямой ход по всем неизвестным уравнениям
    for (int k = 0; k < n; k++) {
        // Исключение x_k из строк k+1...n-1
        double pivot = a[k * n + k];
        for (int i = 0; i < n; i++) {
            // Из урвнения строки i вычитается уравнение k
            double lik = a[i * n + k] / pivot;
            for (int j = 0; j < n; j++)
                a[i * n + j] -= lik * a[k * n + j];
            b[i] -= lik * b[k];
        }
    }

    // Обратный ход
    for (int k = n - 1; k >= 0; k--) {
        x[k] = b[k];
        for (int i = k + 1; i < n; i++)
            x[k] -= a[k * n + i] * x[i];
        x[k] /= a[k * n + k];
    }

}
