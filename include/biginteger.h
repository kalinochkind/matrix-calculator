#ifndef _BIGINTEGER_H
#define _BIGINTEGER_H

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

enum compare_t
{
    CMP_EQUAL = 0, CMP_LESS = -1, CMP_GREATER = 1
};

class invalid_number_error: public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

class zero_division_error: public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

template<class T>
inline const std::string toString(const T &a)
{
    std::ostringstream os;
    os << a;
    return os.str();
}


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

#endif