#ifndef _RATIONAL_H
#define _RATIONAL_H

#include "biginteger.h"

class Rational
{
private:
    BigInteger _numerator, _denominator;

    void normalize();

    void expand(const BigInteger &a);

public:
    Rational(const BigInteger &a): _numerator(a), _denominator(1) {};

    Rational(long long a): _numerator(a), _denominator(1) {};

    Rational(): _numerator(0), _denominator(1) {};

    Rational(const BigInteger &a, const BigInteger &b): _numerator(a), _denominator(b)
    {
        normalize();
    }

    explicit Rational(const std::string &s);

    compare_t compare(const Rational &a) const
    {
        return (_numerator * a._denominator - a._numerator * _denominator)._compare(0);
    }

    const BigInteger numerator() const
    {
        return _numerator;
    }

    void numerator(const BigInteger &a)
    {
        _numerator = a;
        normalize();
    }

    const BigInteger denominator() const
    {
        return _denominator;
    }

    void denominator(const BigInteger &a)
    {
        _denominator = a;
        normalize();
    }

    std::string toString() const;

    explicit operator int() const
    {
        return int(BigInteger(*this));
    }

    explicit operator const BigInteger() const
    {
        if(_denominator == 1)
            return _numerator;
        return _numerator / _denominator;
    }

    const Rational operator+() const
    {
        return *this;
    }

    const Rational operator-() const
    {
        Rational temp(*this);
        temp._numerator = -temp._numerator;
        return temp;
    }

    explicit operator bool() const
    {
        return bool(_numerator);
    }

    Rational &operator+=(const Rational &a);

    Rational &operator-=(const Rational &a)
    {
        return *this += -a;
    }

    Rational &operator*=(const Rational &a);

    Rational &operator/=(const Rational &a);

    const std::string asDecimal(size_t precision) const;

};

inline bool operator<(const Rational &a, const Rational &b)
{
    return a.compare(b) == CMP_LESS;
}

inline bool operator==(const Rational &a, const Rational &b)
{
    return a.compare(b) == CMP_EQUAL;
}

inline bool operator>(const Rational &a, const Rational &b)
{
    return a.compare(b) == CMP_GREATER;
}

inline bool operator<=(const Rational &a, const Rational &b)
{
    return a.compare(b) != CMP_GREATER;
}

inline bool operator>=(const Rational &a, const Rational &b)
{
    return a.compare(b) != CMP_LESS;
}

inline bool operator!=(const Rational &a, const Rational &b)
{
    return a.compare(b) != CMP_EQUAL;
}

const Rational operator+(const Rational &a, const Rational &b);

const Rational operator-(const Rational &a, const Rational &b);

const Rational operator*(const Rational &a, const Rational &b);

const Rational operator/(const Rational &a, const Rational &b);

std::ostream &operator<<(std::ostream &out, const Rational &a);

std::istream &operator>>(std::istream &in, Rational &a);
#endif
