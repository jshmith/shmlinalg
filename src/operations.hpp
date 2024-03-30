#include "matrix.hpp"

/*
Performs the matrix multiplication AB
*/
template <typename T>
Matrix<T> multiply(Matrix<T> &A, Matrix<T> &B);

/*
Performs the matrix addition A+B
*/
template <typename T>
Matrix<T> add(Matrix<T> &A, Matrix<T> &B);

/*
Finds the cofactor matrix of A
*/
template <typename T>
Matrix<T> cofactor(Matrix<T> &A);

/*
Helper for det()
*/
template <typename T>
T detHelper(Matrix<T> &A, T d);

/*
Finds the determinant of matrix A
*/
template <typename T>
T det(Matrix<T> &A);

/*
Finds adj(A)
*/
template <typename T>
Matrix<T> adj(Matrix<T> &A);

/*
Finds the inverse of A
*/
template <typename T>
Matrix<T> inv(Matrix<T> &A);