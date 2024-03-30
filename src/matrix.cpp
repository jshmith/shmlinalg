#include <iostream>
#include <set>
#include "matrix.hpp"
#include "utils.cpp"

// Default constructor
template <typename T>
Matrix<T>::Matrix() : M(0), N(0), data(std::vector<T>()) {
    issquare = false;
}

// Uninitialize constructor
template <typename T>
Matrix<T>::Matrix(int M, int N) : M(M), N(N), data(std::vector<T>(M*N,0)) {
    issquare = (M==N);
}

// Initialize constructor
template <typename T>
Matrix<T>::Matrix(int M, int N, T* data) : M(M), N(N), data(std::vector<T>(data, data + (M*N))) {
    issquare = (M==N);
}

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

// Index method
template <typename T>
T Matrix<T>::Index(size_t row, size_t col) {
    return data[getFlatIndex(M,N,row,col)];
}

// Index overload
template <typename T>
T Matrix<T>::operator[](std::pair<size_t,size_t> index) {
    return data[getFlatIndex(M,N,index.first,index.second)];
}

// Unary minus
template <typename T>
Matrix<T> Matrix<T>::operator-() {
    T outData[M*N];

    for (size_t i = 0; i < M*N; ++i) {
        outData[i] = -data[i];
    }

    Matrix<T> out(M, N, &outData[0U]);
    return out;
}

// Binary minus
template <typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T> B) {
    if ((B.getM() != M) || (B.getN() != N)) {
        throw("Incompatible matrix dimensions.");
    }

    T outData[M*N];

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outData[getFlatIndex(M,N,i,j)] = this->Index(i,j) - B.Index(i,j);
        }
    }

    Matrix<T> out(M, N, outData);
    return out;
}

// Binary division
template <typename T>
Matrix<T> Matrix<T>::operator/(T rhs) {
    T outData[M*N];

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outData[getFlatIndex(M,N,i,j)] = this->Index(i,j) / rhs;
        }
    }

    Matrix<T> out(M, N, outData);
    return out;
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

template <typename T>
bool Matrix<T>::isSquare() {
    return issquare;
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

template <typename T>
void Matrix<T>::printFlat() {
    for (size_t i = 1; i <= M*N; ++i) {
        std::cout << data[i - 1];

        if (i < M*N) {
            std::cout << ", ";
        }
    }
}