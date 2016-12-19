#ifndef _FINITE_H
#define _FINITE_H

#include <iostream>
#include "rational.h"

BigInteger _FINITE_ORDER = 2;

class Finite
{
    BigInteger val;

    const BigInteger binpow(BigInteger a, BigInteger p) const
    {
        BigInteger t = 1;
        while (p)
        {
            if (p.odd())
                t = (t * a) % _FINITE_ORDER;
            a = (a * a) % _FINITE_ORDER;
            p /= 2;
        }
        return t;
    }

public:

    Finite(const BigInteger &n): val(n)
    {
        if (val >= _FINITE_ORDER)
            val %= _FINITE_ORDER;
        if (val < 0)
            val = _FINITE_ORDER - ((0ll - val - 1) % _FINITE_ORDER + 1);
    }

    Finite(int n): Finite(BigInteger(n)) {}

    Finite(): val(0) {};

    explicit operator const BigInteger() const
    {
        return val;
    }

    explicit operator int() const
    {
        return int(val);
    }

    explicit operator bool() const
    {
        return bool(val);
    }

    const Finite operator-() const
    {
        return Finite() - *this;
    }

    const Finite operator+() const
    {
        return *this;
    }

    Finite &operator+=(const Finite &a)
    {
        BigInteger nval = val + a.val;
        if (nval >= _FINITE_ORDER)
        {
            nval -= _FINITE_ORDER;
        }
        val = nval;
        return *this;
    }

    const Finite operator+(const Finite &a) const
    {
        Finite temp(*this);
        return temp += a;
    }

    Finite &operator-=(const Finite &a)
    {
        BigInteger nval = val - a.val;
        if (nval < 0)
        {
            nval += _FINITE_ORDER;
        }
        val = nval;
        return *this;
    }

    const Finite operator-(const Finite &a) const
    {
        Finite temp(*this);
        return temp -= a;
    }

    Finite &operator*=(const Finite &a)
    {
        BigInteger nval = val * a.val;
        if (nval >= _FINITE_ORDER)
        {
            nval %= _FINITE_ORDER;
        }
        val = nval;
        return *this;
    }

    const Finite operator*(const Finite &a) const
    {
        Finite temp(*this);
        return temp *= a;
    }

    const Finite inverse() const
    {
        auto p = ext_gcd(val, _FINITE_ORDER);
        return Finite(p.first);
    }

    Finite operator/=(const Finite &a)
    {
        return *this *= a.inverse();
    }

    const Finite operator/(const Finite &a) const
    {
        return *this * a.inverse();
    }

};


bool operator==(const Finite &a, const Finite &b)
{
    return BigInteger(a) == BigInteger(b);
}

bool operator!=(const Finite &a, const Finite &b)
{
    return BigInteger(a) != BigInteger(b);
}


std::ostream &operator<<(std::ostream &out, const Finite &a)
{
    return out << BigInteger(a);
}

std::istream &operator>>(std::istream &in, Finite &a)
{
    BigInteger t;
    in >> t;
    a = t;
    return in;
}

#endif