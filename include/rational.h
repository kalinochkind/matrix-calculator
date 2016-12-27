#ifndef _RATIONAL_H
#define _RATIONAL_H

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

template<class T>
inline const std::string toString(const T &a)
{
    std::ostringstream os;
    os << a;
    return os.str();
}

enum compare_t
{
    CMP_EQUAL = 0, CMP_LESS = -1, CMP_GREATER = 1
};

class zero_division_error: public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

class BigInteger
{
private:
    std::vector<unsigned> digits;
    bool negative;

    void normalize();

    void shiftRight(size_t a);

    void multiply(unsigned a);

    const BigInteger karatsuba(const BigInteger &a) const;

    compare_t _compareAbs(const BigInteger &a) const;

    void _divide(int a, int &remainder);

    const BigInteger _divide_long(BigInteger a);

public:

    static const unsigned BLOCK_MOD = 1000000000;
    static const size_t BLOCK_SIZE = 9;

    BigInteger(): digits(), negative(false) {};

    BigInteger(long long a);

    explicit BigInteger(const std::string &s);

    compare_t _compare(const BigInteger &a) const;

    explicit operator bool() const
    {
        return !digits.empty();
    }

    explicit operator int() const
    {
        return digits.empty() ? 0 : (digits[0] * (negative ? -1 : 1));
    }

    std::string toString() const;

    const BigInteger operator-() const
    {
        BigInteger a(*this);
        a.negative = !a.negative;
        return a;
    }

    const BigInteger operator+() const
    {
        return *this;
    }

    BigInteger &operator+=(const BigInteger &a);

    BigInteger &operator*=(BigInteger a);

    BigInteger &operator/=(const BigInteger &a);

    BigInteger &operator%=(const BigInteger &a);

    BigInteger &operator-=(const BigInteger &a);

    BigInteger &operator++()
    {
        return *this += 1;
    }

    BigInteger &operator--()
    {
        return *this -= 1;
    }

    const BigInteger operator++(int)
    {
        BigInteger a = *this;
        ++*this;
        return a;
    }

    const BigInteger operator--(int)
    {
        BigInteger a = *this;
        --*this;
        return a;
    }

    bool odd() const
    {
        return digits.size() && digits[0] % 2;
    }

    void swap(BigInteger &a)
    {
        digits.swap(a.digits);
        std::swap(negative, a.negative);
    }
};


const BigInteger operator+(const BigInteger &a, const BigInteger &b);

const BigInteger operator-(const BigInteger &a, const BigInteger &b);
const BigInteger operator*(const BigInteger &a, const BigInteger &b);

const BigInteger operator/(const BigInteger &a, const BigInteger &b);
const BigInteger operator%(const BigInteger &a, const BigInteger &b);
std::ostream &operator<<(std::ostream &out, const BigInteger &a);
std::istream &operator>>(std::istream &in, BigInteger &a);

inline bool operator<(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) == CMP_LESS;
}

inline bool operator==(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) == CMP_EQUAL;
}

inline bool operator>(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) == CMP_GREATER;
}

inline bool operator<=(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) != CMP_GREATER;
}

inline bool operator>=(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) != CMP_LESS;
}

inline bool operator!=(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) != CMP_EQUAL;
}

inline const BigInteger abs(const BigInteger &a)
{
    return (a < 0) ? -a : a;
}

const BigInteger abs(const BigInteger &a);

const BigInteger gcd(BigInteger a, BigInteger b);

std::pair<BigInteger, BigInteger> ext_gcd(BigInteger a, BigInteger b);


class Rational;

const Rational operator+(const Rational &a, const Rational &b);

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
