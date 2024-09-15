#include <iostream>
#include "matrix.hpp"
#include "operations.hpp"
#include <memory>

MatrixPtr<double> makeMatrix() {
    // Define what the data looks like, either as a vector or a C++ array
    // Rows of a Matrix are contiguous in memory.
    std::vector<double> data({1, 2, 3, 1, 2, 3, 1, 2, 3}); // Underlying data
    size_t M = 3; // Number of rows
    size_t N = 3; // Number of columns
    return Matrix<double>::create(M, N, data);
}

MatrixPtr<double> makeVector() {
    std::vector<double> data({-1, 0, 1}); // Underlying data
    size_t M = 3; // Number of rows
    size_t N = 1; // Number of columns
    return Matrix<double>::create(M, N, data);
}

MatrixPtr<double> makeHermitianMatrix() {
    std::vector<double> data({2, -1, -1, 2});
    size_t M = 2;
    return Matrix<double>::create(M, M, data);
}

int main() {
    /*
    Example 1: Create and print a matrix
    */
    std::cout << "****** Example 1 ******" << std::endl;
    MatrixPtr<double> A = makeMatrix();
    A->print();
    std::cout << std::endl << std::endl;

    /*
    Example 2: Multiply matrices and print the result.
    */
    std::cout << "****** Example 2 ******" << std::endl;
    MatrixPtr<double> B = multiply(A, A);
    B->print();
    std::cout << std::endl << std::endl;

    /*
    Example 3: Make random matrix, invert it, multiply the two.
    */
    std::cout << "****** Example 3 ******" << std::endl;
    size_t M = 3;
    MatrixPtr<double> C = randmat<double>(M, M);
    MatrixPtr<double> Cinv = inv(C);
    MatrixPtr<double> D = multiply(Cinv, C);
    D->print();
    std::cout << std::endl << std::endl;

    /*
    Example 4: Make diagonal matrix from a vector.
    */
    std::cout << "****** Example 4 ******" << std::endl;
    MatrixPtr<double> vec = makeVector();
    vec->print();
    MatrixPtr<double> diagMat = Diagonal(vec);
    std::cout << "-----------------" << std::endl;
    diagMat->print();
    std::cout << std::endl << std::endl;

    /*
    Example 5: QR Decomposition.
    */
    std::cout << "****** Example 5 ******" << std::endl;
    size_t N = 2;
    std::cout << "-------Starting Matrix--------" << std::endl;
    MatrixPtr<double> E = randmat<double>(N, N);
    E->print();

    auto QR = qr(E);
    auto Q = QR.first;
    auto R = QR.second;
    std::cout << "-------Q--------" << std::endl;
    Q->print();
    std::cout << "--------R--------" << std::endl;
    R->print();
    std::cout << "--------Product QR--------" << std::endl;
    auto prod = multiply(Q, R);
    prod->print();
    std::cout << std::endl << std::endl;

    /*
    Example 6: Eigenvalue Decomposition.
    */
    std::cout << "****** Example 6 ******" << std::endl;
    MatrixPtr<double> H = makeHermitianMatrix();
    std::cout << "--------Starting Matrix--------" << std::endl;
    H->print();

    auto EIG = eig(H);
    auto L = EIG.first;
    auto V = EIG.second;

    std::cout << "--------L--------" << std::endl;
    L->print();

    std::cout << "--------V--------" << std::endl;
    V->print();

    std::cout << "--------Product VLV'--------" << std::endl;
    auto prod2 = multiply(V, multiply(Diagonal(L), transpose(V)));
    prod2->print();
    return 0;
}