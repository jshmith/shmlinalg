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

/*
Returns the Matrix A except for the row and col.
*/
template <typename T>
Matrix<T> cofactor(Matrix<T> &A, size_t row, size_t col) {
    if (row >= A.getM() || col >= A.getN()) {
        throw("Row or column is out of bounds for the matrix.");
    }

    if (A.getM() == 1 || A.getN() == 1) {
        throw("Cofactor not well defined for a vector.");
    }

    size_t outM = A.getM()-1;
    size_t outN = A.getN()-1;

    T outData[outM * outN];

    size_t pastCol = 0;
    size_t pastRow = 0;
    for (size_t i = 0; i < outM; ++i) {
        pastCol = 0;
        if (i >= row) {
            pastRow = 1;
        }

        for (size_t j = 0; j < outN; ++j) {
            if (j >= col) {
                pastCol = 1;
            }
            outData[outN*i + j] = A.Index(i+pastRow,j+pastCol);
        }
    }

    Matrix<T> out(outM, outN, outData);

    return out;
}

/*
Helper for det()
*/
template <typename T>
T detHelper(Matrix<T> &A, T d) {
    // Base case
    if (A.getM() == 1 && A.getN() == 1) {
        return A.Index(0,0);
    }
    else {
        T sign = -1;
        // Iterate over first row and find the determinant
        // for the cofactor of each element
        for (size_t j = 0; j < A.getN(); ++j) {
            auto C = cofactor(A, 0, j);
            sign *= -1;
            d += A.Index(0,j)*sign*detHelper(C, 0.0);
        }
        return d;
    }
}

/*
Finds the determinant of matrix A
*/
template <typename T>
T det(Matrix<T> &A) {
    // Validate A is square
    if (!A.isSquare()) {
        throw("Matrix must be square.");
    }

    return detHelper(A, 0.0);
}

/*
Finds adj(A)
*/
template <typename T>
Matrix<T> adj(Matrix<T> &A) {
    // Validate A is square
    if (!A.isSquare()) {
        throw("Matrix must be square.");
    }

    size_t dimSize = A.getM();
    T outData[dimSize*dimSize];

    for (size_t i = 0; i < dimSize; ++i) {
        for (size_t j = 0; j < dimSize; ++j) {
            auto C = cofactor(A, i, j);
            T d = det(C);
            T sign = ((i+j)%2) ? -1 : 1;
            size_t ind = getFlatIndex(dimSize,dimSize,j,i);
            outData[ind] = sign*d;
        }
    }

    Matrix<T> out(dimSize,dimSize,outData);

    return out;
}

/*
Finds the inverse of A
*/
template <typename T>
Matrix<T> inv(Matrix<T> &A) {
    // Validate A is square
    if (!A.isSquare()) {
        throw("Matrix must be square.");
    }

    T d = det(A);

    if (d == 0) {
        throw("Matrix is singular. Inverse is not defined.");
    }

    auto Adj = adj(A);

    return Adj/d;
}