#include "matrix.hpp"
#include "operations.hpp"
#include "utils.hpp"
#include <random>
#include <cmath>
#include <cassert>

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

    return Matrix<T>::create(mA, nB, &outData[0U]);
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
            outData[getFlatIndex(mA,nA,i,j)] = A->Index(i, j) + B->Index(i, j);
        }
    }

    // Matrix<T> C(mA, nA, &outData[0U]);
    return Matrix<T>::create(mA,nA,&outData[0U]);
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
    return Matrix<T>::create(N,M,&outData[0U]);
}

/*
Returns the length-N standard basis vector for the i-th dimension
*/
template <typename T>
MatrixPtr<T> standardBasis(size_t N, size_t i) {
    T outData[N];

    for (size_t j = 0; j < N; ++j) {
        if (j != i) {
            outData[j] = 0;
        }
        else {
            outData[j] = 1;
        }
    }

    return Matrix<T>::create(N,1,&outData[0U]);
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

    return Matrix<T>::create(outN,outM,&outData[0U]);
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

    return Matrix<T>::create(dimSize,dimSize,&outData[0U]);
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

    return Matrix<T>::create(Adj/d);
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
    return Matrix<T>::create(M,M,&outData[0U]);
}

/*
Finds the QR decomposition of A
*/
template <typename T>
std::pair<MatrixPtr<T>, MatrixPtr<T>> qr(MatrixPtr<T> A) {
    size_t M = A->getM();
    size_t N = A->getN();

    auto Q = eye<T>(M);
    MatrixPtr<T> R = Matrix<T>::create(*A);

    // Iterate over columns of A
    for (size_t j = 0; j < N; ++j) {
        // Find the Householder reflection for the column
        T Rjj = R->Index(j,j);
        MatrixPtr<T> col = Matrix<T>::create(R->Index(j, M-1, j, j));
        T normCol = col->norm();
        int s = (Rjj < 0) ? 1 : -1;
        T u1 = Rjj - s*normCol;

        MatrixPtr<T> w = Matrix<T>::create(*col / u1);
        w->setElem(0, 0, 1);
        T tau = -s * u1 / normCol;

        auto Rslice = Matrix<T>::create(R->Index(j, M-1, 0, N-1));
        auto Qslice = Matrix<T>::create(Q->Index(0, M-1, j, N-1));
        auto wT = transpose(w);

        auto outerW = multiply(w, wT);

        // H is tau * w * w'
        auto H = Matrix<T>::create(*(multiply(w, wT)) * tau);

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
    return std::make_pair(Matrix<T>::create(*Q), Matrix<T>::create(*R));
}

/*
Checks if the matrix A is "close enough" to
being upper triangular
*/
template <typename T>
bool isUpperTriangular(MatrixPtr<T> A, T tol) {
    size_t M = A->getM();
    size_t N = A->getN();

    // Validate A is square
    if (!A->isSquare() || (M == 1 && N == 1)) {
        throw("Matrix must be square and non-scalar.");
    }

    // We are interested in elements with index i,j-1
    for (size_t i = 1; i <= M-1; ++i) {
        size_t j = i-1;

        if (std::abs(A->Index(i, j)) > tol) {
            return false;
        }
    }
    return true;
}

/*
Extracts the diagonal of square matrix A
and returns as a row vector
*/
template <typename T>
MatrixPtr<T> diag(MatrixPtr<T> A) {
    // Validate A is square
    if (!A->isSquare()) {
        throw("Matrix must be square.");
    }

    size_t M = A->getM();
    T outData[M];

    for (size_t i = 0; i < M; ++i) {
        outData[i] = A->Index(i, i);
    }

    return Matrix<T>::create(1,M,&outData[0U]);
}

/*
Creates diagonal matrix from C++ array.
*/
template <typename T>
MatrixPtr<T> Diagonal(int M, T* src) {
    MatrixPtr<T> out = Matrix<T>::create(M, M);

    for (size_t i = 0; i < M; ++i) {
        out->setElem(i, i, src[i+(M*i)]);
    }

    return out;
}

/*
Creates diagonal matrix from vector.
*/
template <typename T>
MatrixPtr<T> Diagonal(MatrixPtr<T> A) {
    assert(A->isVector());

    // Initialize output matrix to the correct size.
    size_t matSize = 0;
    if(A->isRowVector()) {
        matSize = A->getN();
    }
    else {
        matSize = A->getM();
    }

    MatrixPtr<T> out = Matrix<T>::create(matSize, matSize);

    for (size_t i = 0; i < matSize; ++i) {
        out->setElem(i, i, A->Index(i));
    }

    return out;
}

/*
Generate a uniform random matrix
*/
template <typename T>
MatrixPtr<T> randmat(size_t M, size_t N) {
    T outData[M*N];

    std::uniform_real_distribution<T> unif(0, 1);
    std::default_random_engine re;

    for (size_t i = 0; i < M*N; ++i) {
        outData[i] = unif(re);
    }

    return Matrix<T>::create(M,N,&outData[0U]);
}

/*
Compute the norm of a matrix or vector
(2-norm or Frobenius)
*/
template <typename T>
T norm(MatrixPtr<T> A) {
    size_t M = A->getM();
    size_t N = A->getN();

    T out = 0;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            out += pow(A->Index(i, j), 2);
        }
    }
    return pow(out, 0.5);
}

/*
Returns a matrix of 1's with the shape M x N
*/
template <typename T>
MatrixPtr<T> ones(size_t M, size_t N) {
    T outData[M*N];

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            outData[i*N + j] = 1;
        }
    }
    return Matrix<T>::create(M,N,&outData[0U]);
}

/*
Normalize A by dividing by its norm
*/
template <typename T>
MatrixPtr<T> normalize(MatrixPtr<T> A) {
    T normA = norm(A);
    return Matrix<T>::create((*A) / normA);
}

/*
Inverse power iterate on A given an initial
eigenvalue guess mu
*/
template <typename T>
MatrixPtr<T> invPowerIter(MatrixPtr<T> A, T mu, T tol, size_t maxIter) {
    // Validate A is square
    if (!A->isSquare()) {
        throw("Matrix must be square.");
    }

    size_t M = A->getM();

    MatrixPtr<T> v = normalize(randmat<T>(M, size_t(1)));
    MatrixPtr<T> vNext = nullptr;

    // Perturb the matrix in the case where mu is 0
    // and the matrix is singular
    Matrix<T> muI = (*(eye<T>(M)) * (mu+1e-6));

    for (size_t i = 0; i < maxIter; ++i) {
        MatrixPtr<T> diff = Matrix<T>::create((*A) - muI);
        MatrixPtr<T> diffInv = inv(diff);

        MatrixPtr<T> prod = multiply(diffInv, v);

        vNext = Matrix<T>::create((*prod) / norm(prod));

        T residual = norm(Matrix<T>::create((*vNext) - (*v)));

        if (residual < tol) {
            break;
        }

        v = vNext;
    }
    return v;
}

/*
Finds the eigen-decomposition of A
*/
template <typename T>
std::pair<MatrixPtr<T>, MatrixPtr<T>> eig(MatrixPtr<T> A) {
    // Validate A is square
    if (!A->isSquare()) {
        throw("Matrix must be square.");
    }

    size_t M = A->getM();

    // Iterate and compute QR until Ak is roughly diagonal
    T tol = 1e-10;
    size_t maxIter = 1000;

    MatrixPtr<T> curA = Matrix<T>::create(*A);
    MatrixPtr<T> Q = nullptr;
    MatrixPtr<T> R = nullptr;
    MatrixPtr<T> eigVals = nullptr;

    // Estimate eigenvalues
    for (size_t i = 0; i < maxIter; ++i) {
        std::pair<std::shared_ptr<Matrix<T>>, std::shared_ptr<Matrix<T>>> QR = qr(curA);
        Q = QR.first;
        R = QR.second;
        curA = multiply(R, Q);
        eigVals = diag(curA);

        if (isUpperTriangular(curA, tol)) {
            curA->print();
            std::cout << "************" << std::endl;
            break;
        }
    }

    MatrixPtr<T> eigVecs = Matrix<T>::create(M, M);

    for (size_t i = 0; i < M; ++i) {
        auto vi = invPowerIter(A, eigVals->Index(0, i), tol, maxIter);
        eigVecs->setElem(0, M-1, i, i, *vi);
    }

    return std::make_pair(eigVals, eigVecs);
}