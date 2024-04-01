#include "matrix.hpp"
#include "operations.hpp"
#include "utils.hpp"

/*
Performs the matrix multiplication AB
*/
template <typename T>
MatrixPtr<T> multiply(MatrixPtr<T> A, MatrixPtr<T> B) {
    // Extract dimensions
    const size_t mA = A->getM();
    const size_t nA = A->getN();

    const size_t mB = B->getM();
    const size_t nB = B->getN();

    // Ensure dimension compatibility
    if (nA != mB) {
        throw("Incompatible matrix dimensions.");
    }

    // Output array data
    size_t nElemOut = mA * nB;
    T outData[nElemOut];

    for (size_t k = 0; k < mA; ++k) {
        for (size_t j = 0; j < nB; ++j) {
            size_t dataIdx = nB*k + j;
            outData[dataIdx] = 0;

            for (size_t i = 0; i < nA; ++i) {
                outData[dataIdx] += (A->Index(k,i)*B->Index(i,j));
            }
        }
    }

    return std::make_shared<Matrix<T>>(mA, nB, &outData[0U]);
}

/*
Performs the matrix addition A+B
*/
template <typename T>
MatrixPtr<T> add(MatrixPtr<T> A, MatrixPtr<T> B) {
    // Extract dimensions
    size_t mA = A->getM();
    size_t nA = A->getN();

    size_t mB = B->getM();
    size_t nB = B->getN();

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
            outData[getFlatIndex(mA,nA,i,j)] = (*A)[loc] + (*B)[loc];
        }
    }

    // Matrix<T> C(mA, nA, &outData[0U]);
    return std::make_shared<Matrix<T>>(mA,nA,outData);
}

/*
Transposes A and returns a new Matrix
*/
template <typename T>
MatrixPtr<T> transpose(MatrixPtr<T> A) {
    size_t M = A->getM();
    size_t N = A->getN();

    size_t numel = M*N;
    T outData[numel];

    // Worry about speed later
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outData[M*j + i] = A->Index(i, j);
        }
    }

    //Matrix<T> out(N,M,outData);
    return std::make_shared<Matrix<T>>(N,M,&outData[0U]);
}

/*
Returns the Matrix A except for the row and col.
*/
template <typename T>
MatrixPtr<T> cofactor(MatrixPtr<T> A, size_t row, size_t col) {
    if (row >= A->getM() || col >= A->getN()) {
        throw("Row or column is out of bounds for the matrix.");
    }

    if (A->getM() == 1 || A->getN() == 1) {
        throw("Cofactor not well defined for a vector.");
    }

    size_t outM = A->getM()-1;
    size_t outN = A->getN()-1;

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
            outData[outN*i + j] = A->Index(i+pastRow,j+pastCol);
        }
    }

    //Matrix<T> out(outM, outN, outData);

    return std::make_shared<Matrix<T>>(outN,outM,outData);
}

/*
Helper for det()
*/
template <typename T>
T detHelper(MatrixPtr<T> A, T d) {
    // Base case
    if (A->getM() == 1 && A->getN() == 1) {
        return A->Index(0,0);
    }
    else {
        T sign = -1;
        // Iterate over first row and find the determinant
        // for the cofactor of each element
        for (size_t j = 0; j < A->getN(); ++j) {
            auto C = cofactor(A, 0, j);
            sign *= -1;
            d += A->Index(0,j)*sign*detHelper(C, 0.0);
        }
        return d;
    }
}

/*
Finds the determinant of matrix A
*/
template <typename T>
T det(MatrixPtr<T> A) {
    // Validate A is square
    if (!A->isSquare()) {
        throw("Matrix must be square.");
    }

    return detHelper(A, 0.0);
}

/*
Finds adj(A)
*/
template <typename T>
MatrixPtr<T> adj(MatrixPtr<T> A) {
    // Validate A is square
    if (!A->isSquare()) {
        throw("Matrix must be square.");
    }

    size_t dimSize = A->getM();
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

    //Matrix<T> out(dimSize,dimSize,outData);

    return std::make_shared<Matrix<T>>(dimSize,dimSize,outData);
}

/*
Finds the inverse of A
*/
template <typename T>
MatrixPtr<T> inv(MatrixPtr<T> A) {
    // Validate A is square
    if (!A->isSquare()) {
        throw("Matrix must be square.");
    }

    T d = det(A);

    if (d == 0) {
        throw("Matrix is singular. Inverse is not defined.");
    }

    auto Adj = *(adj(A));

    return std::make_shared<Matrix<T>>(Adj/d);
}

/*
Returns the identity matrix of size M
*/
template <typename T>
MatrixPtr<T> eye(size_t M) {
    T outData[M*M];

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < M; ++j) {
            T val = 0;
            if (i == j) {
                val = 1;
            }
            outData[i*M + j] = val;
        }
    }
    //Matrix<T> out(M,M,outData);
    return std::make_shared<Matrix<T>>(M,M,&outData[0U]);
}

/*
Finds the QR decomposition of A
*/
template <typename T>
std::pair<MatrixPtr<T>, MatrixPtr<T>> qr(MatrixPtr<T> A) {
    size_t M = A->getM();
    size_t N = A->getN();

    auto Q = eye<T>(M);
    MatrixPtr<T> R = std::make_shared<Matrix<T>>(*A);

    // Iterate over columns of A
    for (size_t j = 0; j < N; ++j) {
        // Find the Householder reflection for the column
        T Rjj = R->Index(j,j);
        MatrixPtr<T> col = std::make_shared<Matrix<T>>(R->Index(j, M-1, j, j));
        T normCol = col->norm();
        int s = (Rjj < 0) ? 1 : -1;
        T u1 = Rjj - s*normCol;

        MatrixPtr<T> w = std::make_shared<Matrix<T>>(*col / u1);
        w->setElem(0, 0, 1);
        T tau = -s * u1 / normCol;

        auto Rslice = std::make_shared<Matrix<T>>(R->Index(j, M-1, 0, N-1));
        auto Qslice = std::make_shared<Matrix<T>>(Q->Index(0, M-1, j, N-1));
        auto wT = transpose(w);

        auto outerW = multiply(w, wT);

        // H is tau * w * w'
        auto H = std::make_shared<Matrix<T>>(*(multiply(w, wT)) * tau);

        // Compute HR and QH for the part of the matrix we care about
        auto HR = multiply(H, Rslice);
        auto QH = multiply(Qslice, H);

        // This is computing R-HR and Q-QH
        // since the formula calls for R(I - tau*w*w')
        // and (I - tau*w*w')Q
        auto diffR = *Rslice - *HR;
        auto diffQ = *Qslice - *QH;

        R->setElem(j, M-1, 0, N-1, diffR);
        Q->setElem(0, M-1, j, N-1, diffQ);
    }
    return std::make_pair(std::make_shared<Matrix<T>>(*Q), std::make_shared<Matrix<T>>(*R));
}