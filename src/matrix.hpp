#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

template <typename T>
class Matrix {
    public:
        Matrix();
        Matrix(int M, int N);
        Matrix(int M, int N, T* src);

        void Transpose();

        T Index(size_t row, size_t col);
        T operator[](std::pair<size_t,size_t> index);
        Matrix operator-();
        Matrix operator-(Matrix B);
        Matrix operator/(T rhs);

        size_t getM();
        size_t getN();
        bool isSquare();

        void print();
        void printFlat();

    private:
        // Size properties
        size_t M;
        size_t N;

        // Booleans
        bool issquare;

        // Data
        std::vector<T> data;
};

#endif // MATRIX_H