#ifndef _MATIX_H
#define _MATIX_H

#include "rational.h"
#include <vector>
#include <iostream>
#include <cassert>

const unsigned strassen_threshold = 32;

class matrix_error: public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

template<class Field>
void _strassen(unsigned N, Field *a, Field *b, Field *res)
{
    if (N <= strassen_threshold)
    {
        for (unsigned i = 0; i < N; ++i)
        {
            for (unsigned j = 0; j < N; ++j)
            {
                for (unsigned q = 0; q < N; ++q)
                {
                    res[i * N + j] += a[i * N + q] * b[q * N + j];
                }
            }
        }
        return;
    }
    Field *buf = new Field[N * N / 2];
    Field *p = new Field[N * N / 4 * 7];
#define AQUAD(x, y) a[(x?N*N/2:0)+(y?N/2:0)+i*N+j]
#define BQUAD(x, y) b[(x?N*N/2:0)+(y?N/2:0)+i*N+j]
#define P(n) p[N*N/4*n+i*N/2+j]
#define FORFOR for (unsigned i = 0; i < N / 2; ++i) for (unsigned j = 0; j < N / 2; ++j)
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 0) + AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 0) + BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(1, 0) + AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 0);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 0);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 1) - BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 2 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(1, 0) - BQUAD(0, 0);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 3 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 0) + AQUAD(0, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 4 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(1, 0) - AQUAD(0, 0);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 0) + BQUAD(0, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 5 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 1) - AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(1, 0) + BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 6 * N * N / 4);
    FORFOR
        {
            res[N * i + j] = P(0) + P(3) - P(4) + P(6);
            res[N * i + N / 2 + j] = P(2) + P(4);
            res[N * N / 2 + N * i + j] = P(1) + P(3);
            res[N * N / 2 + N * i + N / 2 + j] = P(0) - P(1) + P(2) + P(5);
        }
    delete[] buf;
    delete[] p;
#undef AQUAD
#undef BQUAD
#undef P
#undef FORFOR
}

template<class Field>
class Matrix
{
    unsigned M, N;

    template<class Ref>
    class MatrixRow
    {
        Field *arr;
    public:
        MatrixRow(Field *row): arr(row) {};

        Ref operator[](unsigned c)
        {
            return arr[c];
        }
    };

    Field *arr;
    typedef MatrixRow<Field &> row_t;
    typedef MatrixRow<const Field &> crow_t;

    Matrix _multiplyCube(const Matrix &a) const
    {
        assert(N == a.M);
        Matrix res(M, a.N);
        for (unsigned i = 0; i < M; ++i)
        {
            for (unsigned j = 0; j < a.N; ++j)
            {
                for (unsigned q = 0; q < N; ++q)
                {
                    res[i][j] += (*this)[i][q] * a[q][j];
                }
            }
        }
        return res;
    }

    unsigned _gauss(Matrix *ext = nullptr)
    {
        if (ext)
        {
            assert(M == ext->M);
        }
        unsigned row = 0;
        for (unsigned col = 0; col < N; ++col)
        {
            for (unsigned i = row; i < M; ++i)
            {
                if ((*this)[i][col])
                {
                    if (i != row)
                    {
                        for (unsigned j = 0; j < N; ++j)
                        {
                            std::swap((*this)[i][j], (*this)[row][j]);
                            (*this)[i][j] = -(*this)[i][j];
                        }
                        if (ext)
                        {
                            for (unsigned j = 0; j < ext->N; ++j)
                            {
                                std::swap((*ext)[i][j], (*ext)[row][j]);
                                (*ext)[i][j] = -(*ext)[i][j];
                            }
                        }
                    }
                    break;
                }
            }
            if (!(*this)[row][col])
            {
                continue;
            }
            for (unsigned i = row + 1; i < M; ++i)
            {
                Field coeff = (*this)[i][col] / (*this)[row][col];
                for (unsigned j = 0; j < N; ++j)
                {
                    (*this)[i][j] -= (*this)[row][j] * coeff;
                }
                if (ext)
                {
                    for (unsigned j = 0; j < ext->N; ++j)
                        (*ext)[i][j] -= (*ext)[row][j] * coeff;
                }
            }
            ++row;
        }
        return row;
    }


    const Matrix _multiplyStrassen(const Matrix &m) const
    {
        assert(N == m.M);
        Matrix res(M, m.N);
        unsigned mx = std::max(M, std::max(N, m.N));
        while (mx & (mx - 1))
            ++mx;
        Field *a = new Field[mx * mx];
        Field *b = new Field[mx * mx];
        Field *r = new Field[mx * mx];
        for (unsigned i = 0; i < M; ++i)
            for (unsigned j = 0; j < N; ++j)
                a[i * mx + j] = (*this)[i][j];
        for (unsigned i = 0; i < N; ++i)
            for (unsigned j = 0; j < m.N; ++j)
                b[i * mx + j] = m[i][j];
        _strassen(mx, a, b, r);
        for (unsigned i = 0; i < M; ++i)
            for (unsigned j = 0; j < m.N; ++j)
                res[i][j] = r[i * mx + j];
        delete[] a;
        delete[] b;
        delete[] r;
        return res;
    }

public:
    Matrix(): M(0), N(0), arr(nullptr) {};

    Matrix(unsigned height, unsigned width): M(height), N(width), arr(new Field[height * width]) {};

    Matrix(unsigned size): Matrix(size, size) {};

    Matrix(const Matrix &a): Matrix(a.M, a.N)
    {
        for (unsigned i = 0; i < M * N; ++i)
        {
            arr[i] = a.arr[i];
        }
    }

    Matrix &operator=(const Matrix &a)
    {
        if (arr == a.arr)
            return *this;
        M = a.M;
        N = a.N;
        delete[] arr;
        arr = new Field[M * N];
        for (unsigned i = 0; i < M * N; ++i)
        {
            arr[i] = a.arr[i];
        }
        return *this;
    }

    ~Matrix()
    {
        delete[] arr;
    }

    static const Matrix identity(unsigned n)
    {
        Matrix res(n);
        for (unsigned i = 0; i < n; ++i)
        {
            res[i][i] = 1;
        }
        return res;
    }

    static const Matrix fromNumber(Field f)
    {
        Matrix a(1, 1);
        a[0][0] = f;
        return a;
    }

    row_t operator[](unsigned r)
    {
        return row_t(this->arr + r * N);
    }

    crow_t operator[](unsigned r) const
    {
        return crow_t(this->arr + r * N);
    }

    Matrix &operator+=(const Matrix &a)
    {
        if (M != a.M || N != a.N)
            throw matrix_error("Trying to add matrices of different size");
        for (unsigned i = 0; i < N * M; ++i)
        {
            arr[i] += a.arr[i];
        }
        return *this;
    }

    const Matrix operator+(const Matrix &a) const
    {
        Matrix temp(*this);
        return temp += a;
    }

    const Matrix operator-() const
    {
        Matrix temp(*this);
        for (unsigned i = 0; i < N * M; ++i)
        {
            temp.arr[i] = -temp.arr[i];
        }
        return temp;
    }

    const Matrix operator+() const
    {
        return *this;
    }

    Matrix &operator-=(const Matrix &a)
    {
        if (M != a.M || N != a.N)
            throw matrix_error("Trying to subtract matrices of different size");
        for (unsigned i = 0; i < N * M; ++i)
        {
            arr[i] -= a.arr[i];
        }
        return *this;
    }

    const Matrix operator-(const Matrix &a) const
    {
        Matrix temp(*this);
        return temp -= a;
    }

    Matrix &operator*=(Field a)
    {
        for (unsigned i = 0; i < N * M; ++i)
        {
            arr[i] *= a;
        }
        return *this;
    }

    const Matrix operator*(Field a) const
    {
        Matrix temp(*this);
        return temp *= a;
    }

    Matrix &operator*=(const Matrix &a)
    {
        return *this = *this * a;
    }

    const Matrix operator*(const Matrix &a) const
    {
        if (N != a.M)
            throw matrix_error("Trying to multiply matrices of different size");
        if (std::max(M, a.N) < strassen_threshold)
        {
            return _multiplyCube(a);
        }
        else
        {
            return _multiplyStrassen(a);
        }
    }

    const Field det() const
    {
        if (N != M)
            throw matrix_error("Trying to calculate determinant of a non-square matrix");
        Matrix tmp(*this);
        if (tmp._gauss() != N)
            return Field();
        Field ans = tmp[0][0];
        for (unsigned j = 1; j < N; ++j)
        {
            ans *= tmp[j][j];
        }
        return ans;
    }

    const Matrix transposed() const
    {
        Matrix res(N, M);
        for (unsigned i = 0; i < M; ++i)
            for (unsigned j = 0; j < N; ++j)
            {
                res[j][i] = (*this)[i][j];
            }
        return res;
    };

    const Field trace() const
    {
        if (N != M)
            throw matrix_error("Trying to calculate trace of a non-square matrix");
        Field ans;
        for (unsigned i = 0; i < N * N; i += N + 1)
        {
            ans += arr[i];
        }
        return ans;
    }

    const Matrix inverted() const
    {
        Matrix tmp(*this);
        Matrix ext = Matrix::identity(N);
        tmp.inverseExt(ext);
        return ext;
    }

    void inverse()
    {
        *this = inverted();
    }

    unsigned rank() const
    {
        Matrix tmp(*this);
        unsigned a = tmp._gauss();
        return a;
    }

    const std::vector<Field> getRow(unsigned n) const
    {
        std::vector<Field> ans(N);
        for (unsigned i = 0; i < N; i++)
        {
            ans[i] = arr[n * N + i];
        }
        return ans;
    }

    const std::vector<Field> getColumn(unsigned n) const
    {
        std::vector<Field> ans(M);
        for (unsigned i = 0; i < M; i++)
        {
            ans[i] = arr[i * N + n];
        }
        return ans;
    }

    bool operator==(const Matrix &a) const
    {
        if (M != a.M || N != a.N)
            return false;
        for (unsigned i = 0; i < N * M; ++i)
        {
            if (arr[i] != a.arr[i])
                return false;
        }
        return true;
    }

    bool operator!=(const Matrix &a) const
    {
        return !(*this == a);
    }

    unsigned height() const
    {
        return M;
    }

    unsigned width() const
    {
        return N;
    }

    const Matrix submatrix(unsigned x1, unsigned y1, unsigned x2, unsigned y2)
    {
        if (x1 > x2)
            std::swap(x1, x2);
        if (y1 > y2)
            std::swap(y1, y2);
        Matrix res(x2 - x1 + 1, y2 - y1 + 1);
        for (unsigned i = x1; i <= x2; ++i)
        {
            for (unsigned j = y1; j <= y2; ++j)
            {
                res[i - x1][j - y1] = (*this)[i][j];
            }
        }
        return res;
    }

    void inverseExt(Matrix &ext)
    {
        if (N != M)
            throw matrix_error("Trying to inverse a non-square matrix");
        if (M != ext.M)
        {
            throw matrix_error("Invalid use of inverseExt");
        }
        if (_gauss(&ext) != N)
        {
            throw matrix_error("Cannot inverse a singular matrix");
        }
        for (unsigned i = N - 1; i != unsigned(-1); --i)
        {
            Field coeff = (*this)[i][i];
            for (unsigned j = 0; j < N; ++j)
            {
                (*this)[i][j] /= coeff;
            }
            for (unsigned j = 0; j < ext.N; ++j)
            {
                ext[i][j] /= coeff;
            }
            for (unsigned j = i + 1; j < N; j++)
            {
                if (!(*this)[i][j])
                    continue;
                for (unsigned q = 0; q < ext.N; q++)
                {
                    ext[i][q] -= ext[j][q] * (*this)[i][j];
                }
                (*this)[i][j] = Field();
            }
        }
    }

    const Matrix power(int pow)
    {
        if (M != N)
            throw matrix_error("Power is only defined for square matrices");
        Matrix t = identity(width());
        Matrix a(*this);
        unsigned p = abs(pow);
        while (p)
        {
            if (p & 1)
                t *= a;
            p >>= 1;
            a *= a;
        }
        if (pow < 0)
            t.inverse();
        return t;
    }
};

template<class Field>
std::ostream &operator<<(std::ostream &out, const Matrix<Field> &a)
{
    if(!a.width() || !a.height())
        return out;
    for (unsigned i = 0; i < a.height(); ++i)
    {
        for (unsigned j = 0; j < a.width(); ++j)
        {
            out << a[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
};

typedef Matrix<Rational> NumMatrix;

template<class Field>
std::istream &operator>>(std::istream &in, Matrix<Field> &a)
{
    for (unsigned i = 0; i < a.height(); ++i)
    {
        for (unsigned j = 0; j < a.width(); ++j)
        {
            in >> a[i][j];
        }
    }
    return in;
};

#endif
