#ifndef MATHFUNC_H
#define MATHFUNC_H

#include <vector>

void dot(const double *x, const double *y, int size, double &res);

void spMV(const int *ia, const int *ja, const double *a, const double *b, int size, double *res);

void axpy(double a, const double *x, const double *y, int size, double *res);

void scan(const int* input, int* output, int n);

/**
 * @brief Copy one array data to another array
 *
 * @tparam T data type
 * @param y where to copy
 * @param x where to copy from
 * @param size array size
 */
template <typename T>
void arrCopy(T *y, const T *x, const unsigned long size) {
    #pragma omp parallel for default(none) shared(x, y, size)
    for (unsigned int i = 0; i < size; i++) {
        y[i] = x[i];
    }
}


/**
 * @brief Initialize array with some value
 *
 * @tparam T data type
 * @param x array to initialize
 * @param a initialization value
 * @param size array size
 */
template <typename T>
void arrInit(T *x, const T a, const unsigned size) {
    #pragma omp parallel for default(none) shared(x, a, size)
    for (unsigned int i = 0; i < size; i++) {
        x[i] = a;
    }
}

#endif //MATHFUNC_H
