#ifndef _POLYNOM_H
#define _POLYNOM_H

#include "matrix.h"

template<class Field>
class Polynom
{
    Matrix<Field> m;

    void extend(unsigned sz);

    void strip();

    void normalize();

    void multiplyX(unsigned pow);

    compare_t ordinalCmp(const Polynom &a) const;

public:

    Polynom(): m(1, 1) {};

    explicit Polynom(unsigned sz): m(1, sz) {};

    explicit Polynom(const Matrix<Field> &matr): m(matr)
    {
        if(matr.height() != 1)
            throw matrix_error("Polynom: matrix 1*n required");
        strip();
    }

    explicit Polynom(const Field &f): m(1, 1)
    {
        m[0][0] = f;
    }

    explicit Polynom(const std::string &s): m(Matrix<Field>::fromRow(s))
    {
        strip();
    }

    const Matrix<Field> toMatrix(unsigned min_width = 0) const;

    inline void swap(Polynom &a)
    {
        m.swap(a.m);
    }

    int degree() const;

    inline explicit operator bool() const
    {
        return degree() >= 0;
    }

    const Polynom gcd(const Polynom &p) const;

    inline const Polynom operator-() const
    {
        return Polynom(-m);
    }

    inline Polynom &operator*=(Field f)
    {
        m *= f;
        return *this;
    }

    inline const Polynom operator*(Field f) const
    {
        Polynom temp(*this);
        return temp *= f;
    }

    Polynom &operator+=(const Polynom &a);

    inline const Polynom operator+(const Polynom &a) const
    {
        Polynom temp(*this);
        return temp += a;
    }

    inline Polynom &operator-=(const Polynom &a)
    {
        return *this += -a;
    }

    inline const Polynom operator-(const Polynom &a) const
    {
        Polynom temp(*this);
        return temp += -a;
    }

    inline Polynom &operator*=(const Polynom &a)
    {
        return *this = *this * a;
    }

    const Polynom operator*(const Polynom &a) const;

    inline Polynom &operator/=(const Polynom &a)
    {
        return *this = *this / a;
    }

    const Polynom operator/(const Polynom &a) const;

    inline Polynom &operator%=(const Polynom &a)
    {
        return *this = *this % a;
    }

    inline const Polynom operator%(const Polynom &a) const
    {
        Polynom temp = *this / a * a;
        return *this - temp;
    }

    const Polynom diff() const;

    bool isOrdinal() const;

    const Polynom ordinalAdd(const Polynom &a) const;

    const Polynom ordinalMul(const Polynom &a) const;

    const Polynom ordinalSub(const Polynom &a) const;

    const std::pair<Polynom, Polynom> ordinalDiv(const Polynom &a) const;

    const Field valueAt(Field a) const;

    const std::vector<Field> roots() const;

    const Polynom power(const BigInteger &pow) const;

};

#endif
