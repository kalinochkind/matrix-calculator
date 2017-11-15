#ifndef _COMPLEX_H
#define _COMPLEX_H

#include <utility>
#include <string>
#include <cassert>

#include "rational.h"

class Complex
{
    Rational _re, _im;

    const Complex inverse() const;

public:
    Complex(const Rational &a): _re(a), _im(0) {};

    Complex(const BigInteger &a): _re(a), _im(0) {};

    Complex(long long a): _re(a), _im(0) {};

    Complex(): _re(0), _im(0) {};

    Complex(const Rational &a, const Rational &b): _re(a), _im(b) {};

    Complex(const std::string &s);


    const Rational &re() const
    {
        return _re;
    }

    const Rational &im() const
    {
        return _im;
    }

    const Complex operator+() const
    {
        return Complex(_re, _im);
    }

    const Complex operator-() const
    {
        return Complex(-_re, -_im);
    }

    explicit operator bool() const
    {
        return bool(_re) || bool(_im);
    }

    Complex &operator+=(const Complex &a);

    Complex &operator-=(const Complex &a);

    Complex &operator*=(const Complex &a);

    Complex &operator/=(const Complex &a);

    std::string toString() const;

    bool operator<(const Complex &a) const
    {
        return std::make_pair(_re, _im) < std::make_pair(a._re, a._im);
    }

    bool operator>(const Complex &a) const
    {
        return a < *this;
    }

    explicit operator const BigInteger() const
    {
        return BigInteger(_re);
    }

    explicit operator int() const
    {
        return int(_re);
    }

    const BigInteger denominator() const
    {
        BigInteger a(_re.denominator()), b(_im.denominator());
        return a * b / gcd(a, b);
    }
};

inline bool operator==(const Complex &a, const Complex &b)
{
    return a.re() == b.re() && a.im() == b.im();
}

inline bool operator!=(const Complex &a, const Complex &b)
{
    return !(a == b);
}

const Complex operator+(const Complex &a, const Complex &b);

const Complex operator-(const Complex &a, const Complex &b);

const Complex operator*(const Complex &a, const Complex &b);

const Complex operator/(const Complex &a, const Complex &b);

std::ostream &operator<<(std::ostream &out, const Complex &a);

std::istream &operator>>(std::istream &in, Complex &a);

#endif