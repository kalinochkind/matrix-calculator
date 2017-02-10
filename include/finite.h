#ifndef _FINITE_H
#define _FINITE_H

#include <iostream>
#include "biginteger.h"

BigInteger _FINITE_ORDER = 2;

class Finite
{
    BigInteger val;

public:

    Finite(const BigInteger &n): val(n)
    {
        if(val >= _FINITE_ORDER)
            val %= _FINITE_ORDER;
        if(val < 0)
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
        if(nval >= _FINITE_ORDER)
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
        if(nval < 0)
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
        if(nval >= _FINITE_ORDER)
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
        if(!val)
            throw zero_division_error("Trying to invert zero");
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


inline bool operator==(const Finite &a, const Finite &b)
{
    return BigInteger(a) == BigInteger(b);
}

inline bool operator!=(const Finite &a, const Finite &b)
{
    return BigInteger(a) != BigInteger(b);
}

inline std::ostream &operator<<(std::ostream &out, const Finite &a)
{
    return out << BigInteger(a);
}

inline std::istream &operator>>(std::istream &in, Finite &a)
{
    BigInteger t;
    in >> t;
    a = t;
    return in;
}

#endif