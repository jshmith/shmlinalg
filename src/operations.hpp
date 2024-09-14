#include "matrix.hpp"
#include <memory>

/*
Performs the matrix multiplication AB
*/
template <typename T>
MatrixPtr<T> multiply(MatrixPtr<T> A, MatrixPtr<T> B);

/*
Performs the matrix addition A+B
*/
template <typename T>
MatrixPtr<T> add(MatrixPtr<T> A, MatrixPtr<T> B);

/*
Performs the matrix addition A+B
*/
template <typename T>
MatrixPtr<T> add(MatrixPtr<T> A, MatrixPtr<T> B);

/*
Transposes A and returns a new Matrix
*/
template <typename T>
MatrixPtr<T> transpose(MatrixPtr<T> A);

/*
Returns the length-N standard basis vector for the i-th dimension
*/
template <typename T>
MatrixPtr<T> standardBasis(size_t N, size_t i);

/*
Finds the cofactor matrix of A
*/
template <typename T>
MatrixPtr<T> cofactor(MatrixPtr<T> A);

/*
Helper for det()
*/
template <typename T>
T detHelper(MatrixPtr<T> A, T d);

/*
Finds the determinant of matrix A
*/
template <typename T>
T det(Matrix<T> &A);

/*
Finds adj(A)
*/
template <typename T>
MatrixPtr<T> adj(MatrixPtr<T> A);

/*
Finds the inverse of A
*/
template <typename T>
MatrixPtr<T> inv(MatrixPtr<T> A);

/*
Returns the identity matrix of size M
*/
template <typename T>
MatrixPtr<T> eye(size_t M);

/*
Finds the QR decomposition of A
*/
template <typename T>
std::pair<MatrixPtr<T>, MatrixPtr<T>> qr(MatrixPtr<T> A);

/*
Checks if the matrix A is "close enough" to
being upper triangular
*/
template <typename T>
bool isUpperTriangular(MatrixPtr<T> A);

/*
Extracts the diagonal of square matrix A
and returns as a row vector
*/
template <typename T>
MatrixPtr<T> diag(MatrixPtr<T> A);

/*
Generate a random matrix
*/
template <typename T>
MatrixPtr<T> randmat(size_t M, size_t N);

/*
Compute the norm of a matrix or vector
(2-norm or Frobenius)
*/
template <typename T>
T norm(MatrixPtr<T> A);

/*
Returns a matrix of 1's with the shape M x N
*/
template <typename T>
MatrixPtr<T> ones(size_t M, size_t N);

/*
Normalize A by dividing by its norm
*/
template <typename T>
MatrixPtr<T> normalize(MatrixPtr<T> A);

/*
Inverse power iterate on A given an initial
eigenvalue guess mu
*/
template <typename T>
MatrixPtr<T> invPowerIter(MatrixPtr<T> A, T mu, T tol, size_t maxIter);

/*
Finds the eigen-decomposition of A
*/
template <typename T>
std::pair<MatrixPtr<T>, MatrixPtr<T>> eig(MatrixPtr<T> A);