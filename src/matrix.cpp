#include <iostream>
#include <set>
#include "matrix.hpp"
#include "utils.cpp"
#include <cassert>
#include <math.h>

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

// Copy constructor
template <typename T>
Matrix<T>::Matrix(const Matrix<T> &A) {
    M = A.M;
    N = A.N;
    data = A.data;
    issquare = A.issquare;
}

// Create a shared pointer to a Matrix
template <typename T>
template <typename... Args>
MatrixPtr<T> Matrix<T>::create(Args&& ... args) {
    return std::make_shared<Matrix<T>>(args...);
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

// Returns the 2 norm of a vector or the
// Frobenius norm of a matrix.
template <typename T>
T Matrix<T>::norm() {
    size_t numel = M*N;
    T out = 0;
    for (size_t i = 0; i < numel; ++i) {
        out += pow(data[i], 2);
    }
    return pow(out, 0.5);
}

// Index methods
template <typename T>
T Matrix<T>::Index(size_t row, size_t col) {
    if ((row < 0 || row >= M) || (col < 0 || col >= N)) {
        throw("Out of bounds.");
    }
    return data[getFlatIndex(M,N,row,col)];
}

template <typename T>
T Matrix<T>::Index(size_t ind) {
    if (ind < 0 || ind >= (M*N)) {
        throw("Out of bounds.");
    }
    return data[ind];
}

template <typename T>
Matrix<T> Matrix<T>::Index(size_t rowStart, size_t rowStop, size_t colStart, size_t colStop) {
    assert((rowStart <= rowStop) && (colStart <= colStop));
    if ((rowStart < 0 || rowStop >= M) || (colStart < 0 || colStop >= N)) {
        throw("Out of bounds.");
    }

    // Compute output sizes
    size_t outM = rowStop - rowStart + 1;
    size_t outN = colStop - colStart + 1;

    T outData[outM*outN];
    size_t outIdx = 0;
    for (size_t i = rowStart; i <= rowStop; ++i) {
        for (size_t j = colStart; j <= colStop; ++j) {
            outData[outIdx] = data[i*N + j];
            outIdx++;
        }
    }

    Matrix<T> out(outM,outN,outData);
    return out;
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

// Scalar division
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

// Scalar multiplication
template <typename T>
Matrix<T> Matrix<T>::operator*(T rhs) {
    T outData[M*N];

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outData[getFlatIndex(M,N,i,j)] = this->Index(i,j) * rhs;
        }
    }

    Matrix<T> out(M, N, outData);
    return out;
}

// Scalar subtraction
template <typename T>
Matrix<T> Matrix<T>::operator-(T rhs) {
    T outData[M*N];

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outData[getFlatIndex(M,N,i,j)] = this->Index(i,j) - rhs;
        }
    }

    Matrix<T> out(M, N, outData);
    return out;
}

// Scalar addition
template <typename T>
Matrix<T> Matrix<T>::operator+(T rhs) {
    T outData[M*N];

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outData[getFlatIndex(M,N,i,j)] = this->Index(i,j) + rhs;
        }
    }

    Matrix<T> out(M, N, outData);
    return out;
}

// Set functions
template <typename T>
void Matrix<T>::setElem(size_t row, size_t col, T val) {
    if ((row < 0 || row >= M) || (col < 0 || col >= N)) {
        throw("Out of bounds.");
    }
    data[row*N + col] = val;
}

template <typename T>
void Matrix<T>::setElem(size_t rowStart, size_t rowStop, size_t colStart, size_t colStop, Matrix<T> &A) {
    assert((rowStart <= rowStop) && (colStart <= colStop));
    assert((rowStop - rowStart + 1 == A.getM()) && (colStop - colStart + 1 == A.getN()));
    if ((rowStart < 0 || rowStop >= M) || (colStart < 0 || colStop >= N)) {
        throw("Out of bounds.");
    }

    // Compute output sizes
    // size_t outM = rowStop - rowStart + 1;
    // size_t outN = colStop - colStart + 1;

    // T outData[outM*outN];
    size_t aN = A.getN();
    for (size_t i = rowStart; i <= rowStop; ++i) {
        for (size_t j = colStart; j <= colStop; ++j) {
            data[i*N + j] = A.Index(i-rowStart, j-colStart);
        }
    }
}

// Get functions
template <typename T>
size_t Matrix<T>::getM() const {
    return M;
}

template <typename T>
size_t Matrix<T>::getN() const {
    return N;
}

template <typename T>
std::vector<T> Matrix<T>::getData() const {
    return data;
}

template <typename T>
bool Matrix<T>::isSquare() const {
    return issquare;
}

template <typename T>
bool Matrix<T>::isColumnVector() const {
    return (N == 1);
}

template <typename T>
bool Matrix<T>::isRowVector() const {
    return (M == 1);
}

template <typename T>
bool Matrix<T>::isVector() const {
    return isColumnVector() || isRowVector();
}

template <typename T>
Matrix<T> Matrix<T>::col(size_t j) {
    if (j >= N) {
        throw("Invalid column for this matrix size.");
    }

    size_t startIdx = j;
    size_t inc = N;

    T outData[M];

    for (size_t i = 0; i < M; ++i) {
        size_t idx = startIdx + inc*i;
        outData[i] = data[idx];
    }

    Matrix<T> out = Matrix<T>(M,1,outData);
    return out;
}

template <typename T>
Matrix<T> Matrix<T>::row(size_t i) {
    if (i >= M) {
        throw("Invalid row for this matrix size.");
    }

    size_t startIdx = i*N;
    size_t inc = 1;

    T outData[N];

    for (size_t i = 0; i < N; ++i) {
        size_t idx = startIdx + inc*i;
        outData[i] = data[idx];
    }

    Matrix<T> out = Matrix<T>(1,M,outData);
    return out;
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