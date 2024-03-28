#include <vector>

template <typename T>
class Matrix {
    public:
        Matrix();
        Matrix(int M, int N);
        Matrix(int M, int N, T* src);

        void Transpose();

        size_t getM();
        size_t getN();

        void print();

    private:
        // Size properties
        size_t M;
        size_t N;

        // Data
        std::vector<T> data;
};