#pragma once

#include "matrix.hpp"
#include "matrix_utilities.hpp"
#include "utils.hpp"

#include <vector>
#include <memory>
#include <iostream>
#include <set>
#include <cassert>
#include <math.h>

struct ForwardPassOutput {
    MatrixPtr<double> output; /* (N, Dout) */
    MatrixPtr<double> X; /* (N, Din) */
    MatrixPtr<double> W; /* (Din, Dout) */
    MatrixPtr<double> b; /* (Dout,) */
};

struct BackwardPassOutput {
    MatrixPtr<double> dX; /* (N, Din) */
    MatrixPtr<double> dW; /* (Din, Dout) */
    MatrixPtr<double> db; /* (Dout,) */
};

// Abstract layer class
class Layer {
    public:
        Layer(MatrixPtr<double> X, size_t dIn, size_t dOut)
        : X(X), W(randn(dIn, dOut)), b(randn(1, dOut)) {}

        virtual ForwardPassOutput& forward();
        
        virtual BackwardPassOutput& backward();
    
    protected:
        MatrixPtr<double> X;
        MatrixPtr<double> W;
        MatrixPtr<double> b;
};

// Concrete, fully-connected layer class
class FCLayer : Layer {
    public:
        FCLayer(MatrixPtr<double> X, size_t dIn, size_t dOut) : Layer(X, dIn, dOut) {}

        ForwardPassOutput& forward() {
            
        }

        BackwardPassOutput& backward() {

        }
};

// Concrete, ReLU layer class
class ReLULayer : Layer {
    public:
        ReLULayer(MatrixPtr<double> X) : Layer(X, 0, 0) {}

        ForwardPassOutput& forward() {
            /* TODO: Implement an binary, elementwise maximum operator for matrices */
        }

        BackwardPassOutput& backward() {
            /* TODO: Implement comparison operators */
        }
};

// Fully connected layer, forward pass
MatrixPtr<double> fc_forward(MatrixPtr<double> X /* Input data (N, Din) */,
                             MatrixPtr<double> w /* Weights (Din, Dout) */,
                             MatrixPtr<double> b /* Biases (Dout,) */);

// Fully connected layer, backward pass
MatrixPtr<double> fc_forward(MatrixPtr<double> diff /* Derivative (N, Dout) */,
                             MatrixPtr<double> X /* Input data (N, Din) */,
                             MatrixPtr<double> w /* Weights (Din, Dout) */,
                             MatrixPtr<double> b /* Biases (Dout,) */);

// ReLU layer, forward pass
MatrixPtr<double> relu_forward(MatrixPtr<double> X /* Input, any shape */);

// ReLU layer, backward pass
MatrixPtr<double> relu_backward(MatrixPtr<double> X /* Input, any shape */);

