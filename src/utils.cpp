#include "utils.hpp"

size_t getFlatIndex(size_t M, size_t N, size_t row, size_t col) {
    (void)M;
    return row*N + col;
}