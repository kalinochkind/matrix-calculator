#ifndef _FINITE_H
#define _FINITE_H

#include <iostream>

int _FINITE_ORDER = 2;

class Finite
{
    int val;

    long long binpow(long long a, int p) const
    {
        long long t = 1;
        while (p)
        {
            if (p & 1)
                t = (t * a) % _FINITE_ORDER;
            a = (a * a) % _FINITE_ORDER;
            p >>= 1;
        }
        return t;
    }

public:

    Finite(int n): val(n)
    {
        if (val >= _FINITE_ORDER)
            val %= _FINITE_ORDER;
        if (val < 0)
            val = _FINITE_ORDER - ((0ll - val - 1) % _FINITE_ORDER + 1);
    };

    Finite(): val(0) {};

    explicit operator int() const
    {
        return val;
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
        long long nval = val + 0ll + a.val;
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
        long long nval = val + 0ll - a.val;
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
        long long nval = val * 1ll * a.val;
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
        return Finite(binpow(val, _FINITE_ORDER - 2));
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
    return int(a) == int(b);
}

bool operator!=(const Finite &a, const Finite &b)
{
    return int(a) != int(b);
}


std::ostream &operator<<(std::ostream &out, const Finite &a)
{
    return out << int(a);
}

std::istream &operator>>(std::istream &in, Finite &a)
{
    int t;
    in >> t;
    a = t;
    return in;
}

#endif