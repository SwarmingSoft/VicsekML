//v11 23.06.2015

#include "MyMatrix.hpp"

RealMatrix::RealMatrix(void) : buffer()
{

}

RealMatrix::RealMatrix(unsigned int rows, unsigned int columns) : buffer(columns, rows)
{

}

unsigned int RealMatrix::GetRows() const
{
    return buffer.GetHeight();
}

unsigned int RealMatrix::GetColumns() const
{
    return buffer.GetWidth();
}

VectorMatrix RealMatrix::operator*(const VectorMatrix &rhs) const
{
    // *this * rhs
    assert(GetColumns() == rhs.GetRows());

    VectorMatrix ret(GetRows(), rhs.GetColumns());
    unsigned int r, c, l;
    Vector value;
    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < rhs.GetColumns(); ++c)
        {
            value.x = real(0); value.y = real(0);
            for (l = 0; l < GetColumns(); ++l)
            {
                value += (*this)(r, l) * rhs(l, c);
                //std::tcout << Get(k, j) << "*" << rhs.Get(i, k) << "+";
            }
            //std::tcout << "=" << value << " to " << i << "," << j;
            ret(r, c) = value;
            //std::tcout << "=> " << ret.Get(i, j) << std::endl;
        }
    }
    return ret;
}

RealMatrix RealMatrix::operator*(const RealMatrix &rhs) const
{
    // *this * rhs
    assert(GetColumns() == rhs.GetRows());

    RealMatrix ret(GetRows(), rhs.GetColumns());
    unsigned int r, c, l;
    real value;
    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < rhs.GetColumns(); ++c)
        {
            value = real(0);
            for (l = 0; l < GetColumns(); ++l)
            {
                value += (*this)(r, l) * rhs(l, c);
                //std::tcout << Get(k, j) << "*" << rhs.Get(i, k) << "+";
            }
            //std::tcout << "=" << value << " to " << i << "," << j;
            ret(r, c) = value;
            //std::tcout << "=> " << ret.Get(i, j) << std::endl;
        }
    }
    return ret;
}

RealMatrix &RealMatrix::operator*=(const RealMatrix &rhs)
{
    *this = *this * rhs;
    return *this;
}

RealMatrix RealMatrix::operator*(real scalar) const
{
    RealMatrix ret(GetRows(), GetColumns());
    unsigned int r, c;

    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < GetColumns(); ++c)
        {
            ret(r, c) = (*this)(r, c) * scalar;
        }
    }
    return ret;
}

RealMatrix RealMatrix::operator-(const RealMatrix &rhs) const
{
    assert(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns());

    RealMatrix ret(GetRows(), GetColumns());
    unsigned int r, c;

    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < GetColumns(); ++c)
        {
            ret(r, c) = (*this)(r, c) - rhs(r, c);
        }
    }

    return ret;
}

RealMatrix RealMatrix::operator+(const RealMatrix &rhs) const
{
    assert(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns());

    RealMatrix ret(GetRows(), GetColumns());
    unsigned int r, c;

    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < GetColumns(); ++c)
        {
            ret(r, c) = (*this)(r, c) + rhs(r, c);
        }
    }
    return ret;
}

RealMatrix &RealMatrix::operator+=(const RealMatrix &rhs)
{
    *this = *this + rhs;
    return *this;
}
RealMatrix &RealMatrix::operator-=(const RealMatrix &rhs)
{
    *this = *this - rhs;
    return *this;
}

RealMatrix::operator real() const
{
    assert(GetRows() == 1 && GetColumns() == 1);
    return (*this)(0, 0);
}

real &RealMatrix::operator()(unsigned int row, unsigned int column)
{
    assert(row < GetRows() && column < GetColumns());
    return *buffer.Get(column, row);
}

const real &RealMatrix::operator()(unsigned int row, unsigned int column) const
{
    assert(row < GetRows() && column < GetColumns());
    return *buffer.Get(column, row);
}

RealMatrix RealMatrix::Identity(unsigned int diagonal)
{
    RealMatrix ret(diagonal, diagonal);
    unsigned int r, c;
    for (r = 0; r < diagonal; ++r)
    {
        for (c = 0; c < diagonal; ++c)
        {
            if (r == c)
            {
                ret(r, c) = real(1);
            }
            else
            {
                ret(r, c) = real(0);
            }
        }
    }
    return ret;
}

RealMatrix RealMatrix::Zero(unsigned int rows, unsigned int columns)
{
    RealMatrix ret(rows, columns);
    unsigned int r, c;
    for (r = 0; r < rows; ++r)
    {
        for (c = 0; c < columns; ++c)
        {
            ret(r, c) = real(0);
        }
    }
    return ret;
}

RealMatrix RealMatrix::Diag(const RealMatrix &m)
{
    assert(m.GetColumns() == 1);

    RealMatrix ret(m.GetRows(), m.GetRows());
    unsigned int r, c;

    for (r = 0; r < ret.GetRows(); ++r)
    {
        for (c = 0; c < ret.GetColumns(); ++c)
        {
            if (r == c)
            {
                ret(r, c) = m(r, 0);
            }
            else
            {
                ret(r, c) = real(0);
            }
        }
    }
    return ret;
}

RealMatrix &RealMatrix::operator=(const RealMatrix &rhs)
{
    buffer = rhs.buffer;
    return *this;
}

RealMatrix &RealMatrix::operator=(const MATRIX &rhs)
{
    assert(buffer.GetWidth() == rhs.GetWidth() && buffer.GetHeight() == rhs.GetHeight());
    unsigned int x, y;

    for (y = 0; y < buffer.GetHeight(); ++y)
    {
        for (x = 0; x < buffer.GetWidth(); ++x)
        {
            *buffer.Get(x, y) = real(*rhs.Get(x, y));
        }
    }
    return *this;
}

bool RealMatrix::operator==(const RealMatrix &rhs) const
{
    return buffer == rhs.buffer;
}

RealMatrix operator*(real scalar, const RealMatrix &rhs)
{
    return rhs*scalar;
}

RealMatrix RealMatrix::Transposed(void) const
{
    RealMatrix ret(GetColumns(), GetRows());
    unsigned int r, c;

    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < GetColumns(); ++c)
        {
            ret(c, r) = (*this)(r, c);
        }
    }

    return ret;
}

RealMatrix &RealMatrix::Transpose(void)
{
    *this = this->Transposed();
    return *this;
}

//###//

VectorMatrix::VectorMatrix(unsigned int rows, unsigned int columns) : buffer(columns, rows)
{

}

unsigned int VectorMatrix::GetRows() const
{
    return buffer.GetHeight();
}

unsigned int VectorMatrix::GetColumns() const
{
    return buffer.GetWidth();
}

RealMatrix VectorMatrix::operator*(const VectorMatrix &rhs) const
{
    // *this * rhs
    assert(GetColumns() == rhs.GetRows());

    RealMatrix ret(GetRows(), rhs.GetColumns());
    unsigned int r, c, l;
    real value;
    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < rhs.GetColumns(); ++c)
        {
            value = real(0);
            for (l = 0; l < GetColumns(); ++l)
            {
                value += (*this)(r, l) * rhs(l, c);
                //std::tcout << Get(k, j) << "*" << rhs.Get(i, k) << "+";
            }
            //std::tcout << "=" << value << " to " << i << "," << j;
            ret(r, c) = value;
            //std::tcout << "=> " << ret.Get(i, j) << std::endl;
        }
    }
    return ret;
}

VectorMatrix VectorMatrix::operator*(const RealMatrix &rhs) const
{
    // *this * rhs
    assert(GetColumns() == rhs.GetRows());

    VectorMatrix ret(GetRows(), rhs.GetColumns());
    unsigned int r, c, l;
    Vector value;
    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < rhs.GetColumns(); ++c)
        {
            value.x = real(0); value.y = real(0);
            for (l = 0; l < GetColumns(); ++l)
            {
                value += (*this)(r, l) * rhs(l, c);
                //std::tcout << Get(k, j) << "*" << rhs.Get(i, k) << "+";
            }
            //std::tcout << "=" << value << " to " << i << "," << j;
            ret(r, c) = value;
            //std::tcout << "=> " << ret.Get(i, j) << std::endl;
        }
    }
    return ret;
}

VectorMatrix VectorMatrix::operator-(const VectorMatrix &rhs) const
{
    assert(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns());

    VectorMatrix ret(GetRows(), GetColumns());
    unsigned int r, c;

    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < GetColumns(); ++c)
        {
            ret(r, c) = (*this)(r, c) - rhs(r, c);
        }
    }
    return ret;
}

Vector &VectorMatrix::operator()(unsigned int row, unsigned int column)
{
    assert(row < GetRows() && column < GetColumns());
    return *buffer.Get(column, row);
}

const Vector &VectorMatrix::operator()(unsigned int row, unsigned int column) const
{
    assert(row < GetRows() && column < GetColumns());
    return *buffer.Get(column, row);
}

VectorMatrix VectorMatrix::Transposed(void) const
{
    VectorMatrix ret(GetColumns(), GetRows());
    unsigned int r, c;

    for (r = 0; r < GetRows(); ++r)
    {
        for (c = 0; c < GetColumns(); ++c)
        {
            ret(c, r) = (*this)(r, c);
        }
    }

    return ret;
}

VectorMatrix &VectorMatrix::Transpose(void)
{
    *this = this->Transposed();
    return *this;
}

bool VectorMatrix::operator==(const VectorMatrix &rhs) const
{
    return buffer == rhs.buffer;
}

RealMatrix MatrixExponential(const RealMatrix &m, unsigned int order)
{
    assert(m.GetRows() == m.GetColumns());
    RealMatrix ret = RealMatrix::Identity(m.GetRows());
    RealMatrix tmp = RealMatrix::Identity(m.GetRows());
    unsigned int denominator = 1;
    for (unsigned int i = 1; i <= order; ++i)
    {
        tmp *= m;
        denominator *= i;
        ret += real(1)/denominator * tmp;
    }

    return ret;
}
