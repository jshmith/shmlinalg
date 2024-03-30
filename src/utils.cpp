#include "utils.hpp"

size_t getFlatIndex(size_t M, size_t N, size_t row, size_t col) {
    return row*N + col;
}