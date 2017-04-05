#ifndef _MATIX_H
#define _MATIX_H

#include <iostream>
#include <cassert>
#include "rational.h"

class matrix_error: public std::runtime_error
{
    using std::runtime_error::runtime_error;
};


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

        inline Ref operator[](unsigned c)
        {
            return arr[c];
        }
    };

    Field *arr;
    typedef MatrixRow<Field &> row_t;
    typedef MatrixRow<const Field &> crow_t;

    static void _strassen(unsigned N, Field *a, Field *b, Field *res);

    Matrix _multiplyCube(const Matrix &a) const;

    const Matrix _multiplyStrassen(const Matrix &m) const;

public:
    Matrix(): M(0), N(0), arr(nullptr) {};

    Matrix(unsigned height, unsigned width): M(height), N(width), arr(new Field[height * width]) {};

    Matrix(unsigned size): Matrix(size, size) {};

    Matrix(const Matrix &a);

    static const Matrix fromRow(const std::string &s);

    Matrix &operator=(const Matrix &a);

    ~Matrix()
    {
        delete[] arr;
    }

    unsigned gauss(Matrix *ext = nullptr);

    static const Matrix identity(unsigned n);

    inline static const Matrix fromNumber(Field f)
    {
        Matrix a(1, 1);
        a[0][0] = f;
        return a;
    }

    inline row_t operator[](unsigned r)
    {
        return row_t(this->arr + r * N);
    }

    inline crow_t operator[](unsigned r) const
    {
        return crow_t(this->arr + r * N);
    }

    Matrix &operator+=(const Matrix &a);

    inline const Matrix operator+(const Matrix &a) const
    {
        Matrix temp(*this);
        return temp += a;
    }

    const Matrix operator-() const;

    inline const Matrix operator+() const
    {
        return *this;
    }

    Matrix &operator-=(const Matrix &a);

    inline const Matrix operator-(const Matrix &a) const
    {
        Matrix temp(*this);
        return temp -= a;
    }

    Matrix &operator*=(Field a);

    inline const Matrix operator*(Field a) const
    {
        Matrix temp(*this);
        return temp *= a;
    }

    inline Matrix &operator/=(Field a)
    {
        return operator*=(Field(1) / a);
    }

    inline const Matrix operator/(Field a) const
    {
        Matrix temp(*this);
        return temp /= a;
    }

    inline Matrix &operator*=(const Matrix &a)
    {
        return *this = *this * a;
    }

    const Matrix operator*(const Matrix &a) const;

    const Field det() const;

    const Matrix transposed() const;

    const Field trace() const;

    const Matrix inverted() const;

    inline unsigned rank() const
    {
        Matrix tmp(*this);
        return tmp.gauss();
    }

    bool operator==(const Matrix &a) const;

    inline bool operator!=(const Matrix &a) const
    {
        return !(*this == a);
    }

    inline unsigned height() const
    {
        return M;
    }

    inline unsigned width() const
    {
        return N;
    }

    const Matrix submatrix(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const;

    void inverseExt(Matrix &ext);

    const Matrix fundamental() const;

    const Matrix partial() const;

    const Matrix power(const BigInteger &pow) const;

    const Matrix joinHorizontal(const Matrix &a) const;

    const Matrix joinVertical(const Matrix &a) const;

    const Matrix charPolynom() const;

    void swap(Matrix &a)
    {
        std::swap(arr, a.arr);
        std::swap(M, a.M);
        std::swap(N, a.N);
    }

    void _write_to_ostream(std::ostream &out) const;

    void _read_from_istream(std::istream &in);
};

template<class Field>
inline std::ostream &operator<<(std::ostream &out, const Matrix<Field> &a)
{
    a._write_to_ostream(out);
    return out;
};

template<class Field>
inline std::istream &operator>>(std::istream &in, Matrix<Field> &a)
{
    a._read_from_istream(in);
    return in;
};

template<class Field>
void makeIntColumns(Matrix<Field> &m);

#endif
