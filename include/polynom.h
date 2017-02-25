#ifndef _POLYNOM_H
#define _POLYNOM_H

#include "matrix.h"

template<class Field>
class Polynom
{
    Matrix<Field> m;

    void extend(unsigned sz)
    {
        if(sz <= m.width())
            return;
        m = Matrix<Field>(1, sz - m.width()).joinHorizontal(m);
    }

    void strip()
    {
        unsigned p = 0;
        while(p < m.width() - 1 && !m[0][p])
            ++p;
        if(p)
            m = m.submatrix(0, p, 0, m.width() - 1);
    }

    void normalize()
    {
        strip();
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

    void multiplyX(unsigned pow)
    {
        extend(m.width() + pow);
        for(unsigned i = 0; i < m.width() - pow; ++i)
            m[0][i] = m[0][i + pow];
        for(unsigned i = m.width() - pow; i < m.width(); ++i)
            m[0][i] = 0;
        strip();
    }

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

    explicit Polynom(const std::string &s): m(1, 1)
    {
        std::istringstream iss;
        iss.str(s);
        Field f;
        std::vector<Field> v;
        while(iss >> f)
            v.push_back(f);
        m = Matrix<Field>(1, v.size());
        for(unsigned i = 0; i < v.size(); ++i)
            m[0][i] = v[i];
        strip();
    }

    const Matrix<Field> toMatrix(unsigned min_width = 0) const
    {
        if(min_width <= m.width())
            return m;
        Polynom t(*this);
        t.extend(min_width);
        return t.m;
    }

    void swap(Polynom &a)
    {
        m.swap(a.m);
    }

    int degree() const
    {
        for(unsigned i = 0; i < m.width(); ++i)
        {
            if(m[0][i])
                return int(m.width() - i) - 1;
        }
        return -1;
    }

    explicit operator bool() const
    {
        return degree() >= 0;
    }

    const Polynom gcd(const Polynom &p) const
    {
        Polynom p1(*this), p2(p), p3;
        while(p2)
        {
            p3 = p1 % p2;
            p3.normalize();
            p1 = p2;
            p2 = p3;
        }
        p1.normalize();
        return p1;
    }

    const Polynom operator-() const
    {
        return Polynom(-m);
    }

    Polynom &operator*=(Field f)
    {
        m *= f;
        return *this;
    }

    const Polynom operator*(Field f) const
    {
        Polynom temp(*this);
        return temp *= f;
    }

    Polynom &operator+=(const Polynom &a)
    {
        extend(a.m.width());
        for(unsigned i = 0; i < a.m.width(); ++i)
        {
            m[0][m.width() - 1 - i] += a.m[0][a.m.width() - 1 - i];
        }
        strip();
        return *this;
    }

    const Polynom operator+(const Polynom &a) const
    {
        Polynom temp(*this);
        return temp += a;
    }

    Polynom &operator-=(const Polynom &a)
    {
        return *this += -a;
    }

    const Polynom operator-(const Polynom &a) const
    {
        Polynom temp(*this);
        return temp += -a;
    }

    Polynom &operator*=(const Polynom &a)
    {
        return *this = *this * a;
    }

    const Polynom operator*(const Polynom &a) const
    {
        Polynom temp(m.width() + a.m.width() - 1);
        for(unsigned i = 0; i < m.width(); ++i)
        {
            for(unsigned j = 0; j < a.m.width(); ++j)
            {
                temp.m[0][temp.m.width() - 1 - i - j] += m[0][m.width() - 1 - i] * a.m[0][a.m.width() - 1 - j];
            }
        }
        temp.strip();
        return temp;
    }

    Polynom &operator/=(const Polynom &a)
    {
        return *this = *this / a;
    }

    const Polynom operator/(const Polynom &a) const
    {
        if(a.m.width() > m.width())
            return Polynom();
        Polynom temp(m.width() - a.m.width() + 1);
        Polynom p1(*this), p2;
        for(unsigned pow = m.width() - a.m.width(); pow != unsigned(-1); --pow)
        {
            p2 = a;
            p2.multiplyX(pow);
            Field f = p1.m[0][0] / p2.m[0][0];
            p1 -= p2 * f;
            temp.m[0][temp.m.width() - pow - 1] = f;
            p1.strip();
            p1.extend(pow + a.m.width() - 1);
        }
        return temp;
    }

    Polynom &operator%=(const Polynom &a)
    {
        return *this = *this % a;
    }

    const Polynom operator%(const Polynom &a) const
    {
        Polynom temp = *this / a * a;
        return *this - temp;
    }

    const Polynom diff() const
    {
        int t = degree();
        if(t < 1)
            return Polynom();
        Polynom ans(t);
        for(int i=1;i<=t;++i)
        {
            ans.m[0][ans.m.width() - i] = m[0][m.width() - i - 1] * i;
        }
        return ans;
    }

};

#endif
