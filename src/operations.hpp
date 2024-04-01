#include "matrix.hpp"
#include <memory>

template <typename T>
using MatrixPtr = std::shared_ptr<Matrix<T>>;

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
Transposes A and returns a new Matrix
*/
template <typename T>
MatrixPtr<T> transpose(MatrixPtr<T> A);

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