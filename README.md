# shmlinalg
This is a C++ Shmlinear Algebra library.

# Usage
The `Matrix` class provides an interface to work with standard matrices. `operations.hpp` is what defines most of the utilities for working iwth `Matrix` objects.

Rather than working directly with `Matrix` instances, the APIs use `MatrixPtr`. These can be made with a call like:
```
MatrixPtr<double> A = Matrix<double>::create(Args... args);
```
See the constructor definitions in `matrix.hpp` for more details on the arguments.

Below is an example on how you might create and print a matrix:
```
size_t M = 2;
auto data = std::vector<double>({1, 2, 3, 4});
auto A = Matrix<double>::create(M, M, data);
A->print();
```

# Try it yourself!
To try this out youself, an `example.cpp` file is provided. To run this, perform the following steps:

1. Run the `make` command in the top-level directory (same level as `example.cpp`)
2. Run `example.exe`

I also encourage you to make modifications to the `example.cpp` file, to try out some things yourself!