#include <iostream>
#include <set>
#include "matrix.hpp"

// Default constructor
template <typename T>
Matrix<T>::Matrix() : M(0), N(0), data(std::vector<T>()) {}

// Uninitialize constructor
template <typename T>
Matrix<T>::Matrix(int M, int N) : M(M), N(N), data(std::vector<T>(M*N,0)) {}

// Initialize constructor
template <typename T>
Matrix<T>::Matrix(int M, int N, T* data) : M(M), N(N), data(std::vector<T>(data, data + (M*N))) {}

// Single matrix operations
template <typename T>
void Matrix<T>::Transpose() {
    // Allocate array for new data
    size_t numel = M*N;
    T newData[numel];

    // Worry about speed later
    for (size_t i = 0; i < numel; ++i) {
        size_t row = i / N;
        size_t col = i % N;
        newData[M*col + row] = data[i];
    }

    size_t tmp = N;
    N = M;
    M = tmp;
    
    data.assign(newData, newData + (M*N));
}

// Index
template <typename T>
T Matrix<T>::Index(size_t row, size_t col) {
    return data[N*row + col];
}

// Get functions
template <typename T>
size_t Matrix<T>::getM() {
    return M;
}

template <typename T>
size_t Matrix<T>::getN() {
    return N;
}

// Print functions
template <typename T>
void Matrix<T>::print() {
    for (size_t i = 1; i <= M*N; ++i) {
        std::cout << data[i - 1];

        if ((i % N) != 0) {
            std::cout << ", ";
        }
        else {
            std::cout << std::endl;
        }
    }
}