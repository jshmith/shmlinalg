#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <memory>

template <typename T>
class Matrix;

template <typename T>
using MatrixPtr = std::shared_ptr<Matrix<T>>;

template <typename T>
class Matrix {
    public:

        Matrix();
        Matrix(int M, int N);
        Matrix(int M, int N, T* src);
        Matrix(const Matrix<T> &A); // Copy constructor
        
        // Static method for creating a shared_ptr to a Matrix
        template <typename... Args>
        static MatrixPtr<T> create(Args&& ... args);

        void Transpose();
        T norm();

        T Index(size_t row, size_t col);
        T Index(size_t ind);
        Matrix Index(size_t rowStart, size_t rowStop, size_t colStart, size_t colStop);
        T operator[](std::pair<size_t,size_t> index);
        Matrix operator-();
        Matrix operator-(Matrix B);
        Matrix operator/(T rhs);
        Matrix operator*(T rhs);
        Matrix operator-(T rhs);
        Matrix operator+(T rhs);

        void setElem(size_t row, size_t col, T val);
        void setElem(size_t rowStart, size_t rowStop, size_t colStart, size_t colStop, Matrix<T> &A);

        size_t getM() const;
        size_t getN() const;
        std::vector<T> getData() const;
        bool isSquare() const;
        bool isColumnVector() const;
        bool isRowVector() const;
        bool isVector() const;

        Matrix col(size_t j);
        Matrix row(size_t i);

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