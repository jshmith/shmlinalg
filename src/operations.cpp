#include "matrix.hpp"
#include "operations.hpp"
#include "utils.hpp"

/*
Performs the matrix multiplication AB
*/
template <typename T>
Matrix<T> multiply(Matrix<T> &A, Matrix<T> &B) {
    // Extract dimensions
    size_t mA = A.getM();
    size_t nA = A.getN();

    size_t mB = B.getM();
    size_t nB = B.getN();

    // Ensure dimension compatibility
    if (nA != mB) {
        throw("Incompatible matrix dimensions.");
    }

    // Output array data
    size_t nElemOut = mA * nB;
    T outData[nElemOut];

    for (size_t k = 0; k < mA; ++k) {
        for (size_t j = 0; j < nB; ++j) {
            size_t dataIdx = mA*k + j;
            outData[dataIdx] = 0;

            for (size_t i = 0; i < nA; ++i) {
                outData[dataIdx] += (A.Index(k,i)*B.Index(i,j));
            }
        }
    }

    Matrix<T> C(mA, nB, &outData[0U]);
    return C;
}

/*
Performs the matrix addition A+B
*/
template <typename T>
Matrix<T> add(Matrix<T> A, Matrix<T> B) {
    // Extract dimensions
    size_t mA = A.getM();
    size_t nA = A.getN();

    size_t mB = B.getM();
    size_t nB = B.getN();

    // Ensure dimension compatibility
    if ((mA != mB) || (nA != nB)) {
        throw("Incompatible matrix dimensions.");
    }

    // Output array data
    size_t nElemOut = mA * nA;
    T outData[nElemOut];

    for (size_t i = 0; i < mA; ++i) {
        for (size_t j = 0; j < nA; ++j) {
            auto loc = std::make_pair(i, j);
            outData[getFlatIndex(mA,nA,i,j)] = A[loc] + B[loc];
        }
    }

    Matrix<T> C(mA, nA, &outData[0U]);
    return C;
}