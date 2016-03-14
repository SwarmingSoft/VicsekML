//v11 23.06.2015

#include <cassert>

#include "constants.hpp"

#include "Buffer.hpp"
#include "Vector2D.hpp"

class VectorMatrix;

class RealMatrix
{
    GenericBuffer_dyn<real, MemoryConst::Dynamic, MemoryConst::Dynamic> buffer;

public:

    RealMatrix(void);
    RealMatrix(unsigned int rows, unsigned int columns);
    unsigned int GetRows() const;
    unsigned int GetColumns() const;
    VectorMatrix operator*(const VectorMatrix &rhs) const;
    RealMatrix operator*(const RealMatrix &rhs) const;
    RealMatrix &operator*=(const RealMatrix &rhs);
    RealMatrix operator*(real scalar) const;
    RealMatrix operator-(const RealMatrix &rhs) const;
    RealMatrix operator+(const RealMatrix &rhs) const;
    RealMatrix &operator+=(const RealMatrix &rhs);
    RealMatrix &operator-=(const RealMatrix &rhs);
    explicit operator real() const;
    real &operator()(unsigned int row, unsigned int column);
    const real &operator()(unsigned int row, unsigned int column) const;
    static RealMatrix Identity(unsigned int diagonal);
    static RealMatrix Zero(unsigned int rows, unsigned int columns);
    static RealMatrix Diag(const RealMatrix &m);
    RealMatrix &operator=(const RealMatrix &rhs);
    RealMatrix &operator=(const MATRIX &rhs);
    bool operator==(const RealMatrix &rhs) const;
    RealMatrix Transposed(void) const;
    RealMatrix &Transpose(void);
};

RealMatrix operator*(real scalar, const RealMatrix &rhs);

class VectorMatrix
{
    GenericBuffer_dyn<Vector, MemoryConst::Dynamic, MemoryConst::Dynamic> buffer;

public:

    VectorMatrix(unsigned int rows, unsigned int columns);
    unsigned int GetRows() const;
    unsigned int GetColumns() const;
    RealMatrix operator*(const VectorMatrix &rhs) const;
    VectorMatrix operator*(const RealMatrix &rhs) const;
    VectorMatrix operator-(const VectorMatrix &rhs) const;
    Vector &operator()(unsigned int row, unsigned int column);
    const Vector &operator()(unsigned int row, unsigned int column) const;
    VectorMatrix Transposed(void) const;
    VectorMatrix &Transpose(void);
    bool operator==(const VectorMatrix &rhs) const;
};

template <typename T> void PrintMatrix(const T &matrix)
{
    unsigned int r, c;

    for (r = 0; r < matrix.GetRows(); ++r)
    {
        for (c = 0; c < matrix.GetColumns(); ++c)
        {
            std::cout << matrix(r, c) << " ";
        }
        std::cout << "\n";
    }
}

RealMatrix MatrixExponential(const RealMatrix &m, unsigned int order);
