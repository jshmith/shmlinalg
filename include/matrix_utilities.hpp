#pragma once

#include "matrix.hpp"

#include <cmath>
#include <random>
#include <vector>

// Generate a standard random Gaussian matrix of size M x N
MatrixPtr<double> randn(size_t M, size_t N) {
    std::normal_distribution d{0.0, 1.0};
    std::random_device rd{};

    // Initialize std::vector of size M * N
    std::vector<double> data(M * N);

    for (size_t i = 0; i < M*N; ++i) {
        data[i] = d(rd);
    }

    return Matrix<double>::create(M, N, data);
}