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

    void multiplyX(unsigned pow)
    {
        extend(m.width() + pow);
        for(unsigned i = 0; i < m.width() - pow; ++i)
            m[0][i] = m[0][i + pow];
        for(unsigned i = m.width() - pow; i < m.width(); ++i)
            m[0][i] = 0;
        strip();
    }

    compare_t ordinalCmp(const Polynom &a) const
    {
        int d = degree(), ad = a.degree();
        if(d < ad)
            return CMP_LESS;
        if(d > ad)
            return CMP_GREATER;
        for(int i=0;i<=d;++i)
        {
            if(m[0][i] < a.m[0][i])
                return CMP_LESS;
            if(m[0][i] > a.m[0][i])
                return CMP_GREATER;
        }
        return CMP_EQUAL;
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

    explicit Polynom(const std::string &s): m(Matrix<Field>::fromRow(s))
    {
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

    bool isOrdinal() const
    {
        for(unsigned i=0;i<m.width();++i)
        {
            if(m[0][i] != BigInteger(m[0][i]))
                return false;
            if(BigInteger(m[0][i]) < 0)
                return false;
        }
        return true;
    }

    const Polynom ordinalAdd(const Polynom &a) const
    {
        if(!isOrdinal() || !a.isOrdinal())
            throw matrix_error("Invalid ordinal");
        if(a.degree() < 0)
            return *this;
        int degdiff = degree() - a.degree();
        if(degdiff < 0)
            return a;
        Polynom x(degree() + 1);
        for(int i=0;i<degdiff;++i)
        {
            x.m[0][i] = m[0][i];
        }
        x.m[0][degdiff] = m[0][degdiff] + a.m[0][0];
        for(unsigned i=degdiff+1;i<x.m.width();++i)
        {
            x.m[0][i] = a.m[0][i - degdiff];
        }
        return x;
    }

    const Polynom ordinalMul(const Polynom &a) const
    {
        if(!isOrdinal() || !a.isOrdinal())
            throw matrix_error("Invalid ordinal");
        int deg = degree();
        if(deg < 0 || a.degree() < 0)
            return Polynom();
        Polynom res(deg + a.degree() + 1);
        for(unsigned i=0;i<a.m.width();++i)
        {
            unsigned pow = a.m.width() - i - 1;
            if(pow > 0)
            {
                Polynom t(a.m[0][i]);
                t.multiplyX(pow + deg);
                res = res.ordinalAdd(t);
            }
            else if(a.m[0][i])
            {
                Polynom t(*this);
                t.m[0][0] *= a.m[0][i];
                res = res.ordinalAdd(t);
            }
        }
        res.strip();
        return res;
    }

    const Polynom ordinalSub(const Polynom &a) const
    {
        if(!isOrdinal() || !a.isOrdinal())
            throw matrix_error("Invalid ordinal");
        int d = degree();
        if(d < 0 || d < a.degree())
            return Polynom();
        Matrix<Field> m(1, d + 1);
        for(int coeff=0;coeff<d+1;++coeff)
        {
            BigInteger l = 0, r = 1;
            for(;;r*=2)
            {
                m[0][coeff] = r - 1;
                Polynom t(m);
                Polynom res = a.ordinalAdd(t);
                compare_t cmp = ordinalCmp(res);
                if(cmp == CMP_EQUAL)
                    return t;
                if(cmp == CMP_LESS)
                    break;
            }
            while(l + 1 < r)
            {
                BigInteger mid = (l + r) / 2;
                m[0][coeff] = mid;
                Polynom t(m);
                Polynom res = a.ordinalAdd(t);
                compare_t cmp = ordinalCmp(res);
                if(cmp == CMP_EQUAL)
                    return t;
                if(cmp == CMP_LESS)
                    r = mid;
                else
                    l = mid;
            }
            m[0][coeff] = l;
        }
        assert(ordinalCmp(a.ordinalAdd(Polynom(m))) == CMP_EQUAL);
        return Polynom(m);
    }

    const std::pair<Polynom, Polynom> ordinalDiv(const Polynom &a) const
    {
        if(!isOrdinal() || !a.isOrdinal())
            throw matrix_error("Invalid ordinal");
        if(a.degree() < 0)
            throw matrix_error("Ordinal division by zero");
        int d = degree();
        if(d < 0 || d < a.degree())
            return {Polynom(), *this};
        Matrix<Field> m(1, d + 1);
        for(int coeff=0;coeff<d+1;++coeff)
        {
            BigInteger l = 0, r = 1;
            for(;;r*=2)
            {
                m[0][coeff] = r - 1;
                Polynom t(m);
                Polynom res = a.ordinalMul(t);
                compare_t cmp = ordinalCmp(res);
                if(cmp == CMP_EQUAL)
                    return {t, Polynom()};
                if(cmp == CMP_LESS)
                    break;
            }
            while(l + 1 < r)
            {
                BigInteger mid = (l + r) / 2;
                m[0][coeff] = mid;
                Polynom t(m);
                Polynom res = a.ordinalMul(t);
                compare_t cmp = ordinalCmp(res);
                if(cmp == CMP_EQUAL)
                    return {t, Polynom()};
                if(cmp == CMP_LESS)
                    r = mid;
                else
                    l = mid;
            }
            m[0][coeff] = l;
        }
        Polynom rem = ordinalSub(a.ordinalMul(Polynom(m)));
        assert(a.ordinalCmp(rem) == CMP_GREATER);
        return {Polynom(m), rem};
    }


};

#endif
