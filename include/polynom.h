#ifndef _POLYNOM_H
#define _POLYNOM_H

#include "matrix.h"

template<class Field>
struct Polynom
{
    Matrix<Field> m;

    void extend(unsigned sz)
    {
        if(sz <= m.width())
            return;
        m = Matrix<Field>(1, sz - m.width()).joinHorizontal(m);
    }

    void normalize()
    {
        unsigned p = 0;
        while(p < m.width() - 1 && !m[0][p])
            ++p;
        if(p)
            m = m.submatrix(0, p, 0, m.width() - 1);
        if(!m[0][0] || m[0][0] == Field(1))
            return;
        m /= m[0][0];
    }

    bool operator<(const Polynom &a) const
    {
        assert(m.width() == a.m.width());
        return degree() < a.degree();
    }

    bool operator>(const Polynom &a) const
    {
        return a < *this;
    }

    bool operator==(const Polynom &a) const
    {
        assert(m.width() == a.m.width());
        return m == a.m;
    }

    bool operator!=(const Polynom &a) const
    {
        return !operator==(a);
    }

    Polynom(): m(1, 1) {};

    Polynom(const Matrix<Field> &matr): m(matr)
    {
        if(matr.height() != 1)
            throw matrix_error("Polynom: matrix 1*n required");
        normalize();
    }

    const Matrix<Field> toMatrix() const
    {
        return m;
    }

    void swap(Polynom &a)
    {
        m.swap(a.m);
    }

    unsigned degree() const
    {
        for(unsigned i = 0; i < m.width(); ++i)
        {
            if(m[0][i])
                return m.width() - i;
        }
        return 0;
    }

    void multiplyX(unsigned pow)
    {
        for(unsigned i = 0; i < m.width() - pow; ++i)
            m[0][i] = m[0][i + pow];
        for(unsigned i = m.width() - pow; i < m.width(); ++i)
            m[0][i] = 0;
    }

    const Polynom gcd(const Polynom &p) const
    {
        Polynom p1(*this), p2(p), p3;
        if(p1.m[0][0] == 0 || p2.m[0][0] == 0)
            return Polynom();
        p1.extend(p2.m.width());
        p2.extend(p1.m.width());
        while(p1 != p2)
        {
            if(p1 < p2)
                p1.swap(p2);
            p3 = p2;
            p3.multiplyX(p1.degree() - p3.degree());
            p1.m -= p3.m;
            p1.normalize();
            p2.normalize();
            p1.extend(p2.m.width());
            p2.extend(p1.m.width());
        }
        return p1;
    }

};

#endif
