#include "polynom.h"
#include "finite.h"

template<class Field>
void Polynom<Field>::extend(unsigned sz)
{
    if(sz > m.width())
        m = Matrix<Field>(1, sz - m.width()).joinHorizontal(m);
}

template<class Field>
void Polynom<Field>::strip()
{
    unsigned p = 0;
    while(p < m.width() - 1 && !m[0][p])
        ++p;
    if(p)
        m = m.submatrix(0, p, 0, m.width() - 1);
}

template<class Field>
void Polynom<Field>::normalize()
{
    strip();
    if(!m[0][0] || m[0][0] == Field(1))
        return;
    m /= m[0][0];
}

template<class Field>
void Polynom<Field>::multiplyX(unsigned pow)
{
    extend(m.width() + pow);
    for(unsigned i = 0; i < m.width() - pow; ++i)
        m[0][i] = m[0][i + pow];
    for(unsigned i = m.width() - pow; i < m.width(); ++i)
        m[0][i] = 0;
    strip();
}

template<class Field>
compare_t Polynom<Field>::ordinalCmp(const Polynom &a) const
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

template<class Field>
const Matrix<Field> Polynom<Field>::toMatrix(unsigned min_width) const
{
    if(min_width <= m.width())
        return m;
    Polynom t(*this);
    t.extend(min_width);
    return t.m;
}

template<class Field>
int Polynom<Field>::degree() const
{
    for(unsigned i = 0; i < m.width(); ++i)
    {
        if(m[0][i])
            return int(m.width() - i) - 1;
    }
    return -1;
}

template<class Field>
const Polynom<Field> Polynom<Field>::gcd(const Polynom &p) const
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

template<class Field>
Polynom<Field> &Polynom<Field>::operator+=(const Polynom &a)
{
    extend(a.m.width());
    for(unsigned i = 0; i < a.m.width(); ++i)
    {
        m[0][m.width() - 1 - i] += a.m[0][a.m.width() - 1 - i];
    }
    strip();
    return *this;
}

template<class Field>
const Polynom<Field> Polynom<Field>::operator*(const Polynom &a) const
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

template<class Field>
const Polynom<Field> Polynom<Field>::operator/(const Polynom &a) const
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

template<class Field>
const Polynom<Field> Polynom<Field>::diff() const
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

template<class Field>
bool Polynom<Field>::isOrdinal() const
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

template<class Field>
const Polynom<Field> Polynom<Field>::ordinalAdd(const Polynom &a) const
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

template<class Field>
const Polynom<Field> Polynom<Field>::ordinalMul(const Polynom &a) const
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

template<class Field>
const Polynom<Field> Polynom<Field>::ordinalSub(const Polynom &a) const
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
    return Polynom(m);
}

template<class Field>
const std::pair<Polynom<Field>, Polynom<Field>> Polynom<Field>::ordinalDiv(const Polynom &a) const
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
    return {Polynom(m), rem};
}

template class Polynom<Rational>;
template class Polynom<Finite>;
