#ifndef _RATIONAL_H
#define _RATIONAL_H

#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <stdexcept>
#include <sstream>

template<class T>
const std::string toString(const T &a)
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

class BigInteger;

const BigInteger operator+(const BigInteger &a, const BigInteger &b);

const BigInteger operator-(const BigInteger &a, const BigInteger &b);

const BigInteger operator*(const BigInteger &a, const BigInteger &b);

const BigInteger operator/(const BigInteger &a, const BigInteger &b);

const BigInteger operator%(const BigInteger &a, const BigInteger &b);

std::ostream &operator<<(std::ostream &out, const BigInteger &a);


class BigInteger
{
private:
    std::vector<unsigned> digits;
    bool negative;

    void normalize()
    {
        while (digits.size() && digits.back() == 0)
            digits.pop_back();
        if (digits.empty())
            negative = false;
    }

    void shiftRight(size_t a)
    {
        digits.resize(digits.size() + a);
        for (size_t i = digits.size() - 1; i >= a; --i)
        {
            digits[i] = digits[i - a];
        }
        for (size_t i = 0; i < a; i++)
        {
            digits[i] = 0;
        }
    }

    void multiply(unsigned a)
    {
        long long carry = 0;
        for (size_t i = 0; carry || i < digits.size(); i++)
        {
            if (i >= digits.size())
                digits.push_back(0);
            long long new_carry = (digits[i] * 1ll * a + carry) / BLOCK_MOD;
            digits[i] = (digits[i] * 1ll * a + carry) % BLOCK_MOD;
            carry = new_carry;
        }
    }

    const BigInteger karatsuba(const BigInteger &a) const
    {
        if (a.digits.size() <= 1)
        {
            BigInteger res(*this);
            res.multiply(a.digits[0]);
            return res;
        }
        size_t maxlen = std::max(digits.size(), a.digits.size());
        maxlen += maxlen % 2;
        BigInteger a1, a2, b1, b2;
        for (size_t i = 0; i < maxlen; i++)
        {
            if (i * 2 < maxlen)
            {
                a2.digits.push_back(i < digits.size() ? digits[i] : 0);
                b2.digits.push_back(i < a.digits.size() ? a.digits[i] : 0);
            }
            else
            {
                a1.digits.push_back(i < digits.size() ? digits[i] : 0);
                b1.digits.push_back(i < a.digits.size() ? a.digits[i] : 0);
            }
        }
        a1.normalize();
        b1.normalize();
        a2.normalize();
        b2.normalize();
        BigInteger a1b1 = a1.karatsuba(b1);
        BigInteger a2b2 = a2.karatsuba(b2);
        BigInteger sum = (a1 + a2).karatsuba(b1 + b2) - a1b1 - a2b2;
        a1b1.shiftRight(maxlen);
        sum.shiftRight(maxlen / 2);
        return a1b1 + sum + a2b2;
    }

    compare_t _compareAbs(const BigInteger &a) const
    {
        if (digits.size() < a.digits.size())
            return CMP_LESS;
        if (digits.size() > a.digits.size())
            return CMP_GREATER;
        for (size_t i = digits.size(); i > 0; --i)
        {
            if (digits[i - 1] < a.digits[i - 1])
                return CMP_LESS;
            if (digits[i - 1] > a.digits[i - 1])
                return CMP_GREATER;
        }
        return CMP_EQUAL;
    }

    void _divide(int a, int &remainder)
    {
        if (!a)
            throw zero_division_error("Integer division by zero");
        if (digits.empty())
        {
            remainder = 0;
            return;
        }
        if (a < 0)
            negative ^= 1;
        long long carry = 0;
        for (auto i = digits.rbegin(); i != digits.rend(); ++i)
        {
            carry = carry * BLOCK_MOD + *i;
            *i = carry / abs(a);
            carry %= abs(a);
        }
        remainder = static_cast<int>(carry);
        normalize();
    }

public:

    static const unsigned BLOCK_MOD = 1000000000;
    static const size_t BLOCK_SIZE = 9;

    BigInteger(): digits(), negative(false) {};

    BigInteger(long long a): digits(), negative(false)
    {
        if (a < 0)
        {
            negative = true;
            a = -a;
        }
        while (a)
        {
            digits.push_back(a % BLOCK_MOD);
            a /= BLOCK_MOD;
        }
    }

    explicit BigInteger(const std::string &s): digits(), negative(false)
    {
        unsigned temp = 0;
        unsigned pow = 1;
        for (auto i = s.crbegin(); i != s.crend(); ++i)
        {
            if (*i == '-')
            {
                negative = true;
                continue;
            }
            temp += pow * (*i - '0');
            pow *= 10;
            if (pow == BLOCK_MOD)
            {
                pow = 1;
                digits.push_back(temp);
                temp = 0;
            }
        }
        digits.push_back(temp);
        normalize();
    }

    compare_t _compare(const BigInteger &a) const
    {
        if (negative && !a.negative)
            return CMP_LESS;
        if (!negative && a.negative)
            return CMP_GREATER;
        return static_cast<compare_t>(_compareAbs(a) * (negative ? -1 : 1));
    }

    explicit operator bool() const
    {
        return !digits.empty();
    }

    explicit operator double() const
    {
        double res = 0;
        for (unsigned i : digits)
        {
            res *= BLOCK_MOD;
            res += i;
        }
        return res;
    }

    explicit operator int() const
    {
        return digits.empty() ? 0 : (digits[0] * (negative ? -1 : 1));
    }

    std::string toString() const
    {
        if (digits.empty())
            return "0";
        std::string res = (negative ? "-" : "");
        for (auto i = digits.rbegin(); i != digits.rend(); ++i)
        {
            std::string temp = ::toString(*i);
            if (i != digits.rbegin())
            {
                while (temp.length() < BLOCK_SIZE)
                {
                    temp = "0" + temp;
                }
            }
            res += temp;
        }
        return res;
    }

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

    BigInteger &operator+=(const BigInteger &a)
    {
        bool sub = negative ^a.negative;
        int rem = 0;
        if (!sub)  // addition
        {
            for (size_t i = 0; i < std::max(a.digits.size(), digits.size()) || rem; ++i)
            {
                if (digits.size() <= i)
                    digits.push_back(0);
                digits[i] = digits[i] + (a.digits.size() > i ? a.digits[i] : 0) + rem;
                rem = digits[i] / BLOCK_MOD;
                digits[i] %= BLOCK_MOD;
            }
        }
        else
        {
            bool rev = _compareAbs(a) < 0;
            const BigInteger &big = (rev ? a : *this), &small = (rev ? *this : a);
            negative = big.negative;
            for (size_t i = 0; i < a.digits.size() || rem; ++i)
            {
                if (digits.size() <= i)
                    digits.push_back(0);
                unsigned dig = (small.digits.size() > i ? small.digits[i] : 0);
                if (big.digits[i] >= dig + rem)
                {
                    digits[i] = big.digits[i] - dig - rem;
                    rem = 0;
                }
                else
                {
                    digits[i] = BLOCK_MOD + big.digits[i] - dig - rem;
                    rem = 1;
                }
            }
        }
        normalize();
        return *this;
    }

    BigInteger &operator*=(BigInteger a)
    {
        bool neg = negative ^a.negative;
        if (_compareAbs(a) == CMP_LESS)
        {
            std::swap(*this, a);
        }
        if (!a.digits.size())
        {
            digits.clear();
        }
        else if (a.digits.size() == 1)
        {
            multiply(a.digits[0]);
        }
        else
        {
            *this = karatsuba(a);
        }
        negative = neg;
        normalize();
        return *this;
    }

    BigInteger &operator/=(const BigInteger &a)
    {
        if (!a)
            throw zero_division_error("Integer division by zero");
        if (a.digits.size() <= 1)
        {
            int temp;
            _divide(a.digits[0], temp);
            return *this;
        }
        BigInteger l(0);
        BigInteger r(*this);
        r.negative = false;
        ++r;
        BigInteger m;
        while ((l + 1)._compareAbs(r) < 0)
        {
            m = (l + r) / 2;
            if (_compareAbs(m * a) >= 0)
            {
                l = m;
            }
            else
            {
                r = m;
            }
        }
        l.negative = a.negative ^ negative;
        l.normalize();
        return *this = l;
    }

    BigInteger &operator%=(const BigInteger &a)
    {
        BigInteger temp(*this);
        temp /= a;
        temp *= a;
        return *this -= temp;
    }

    BigInteger &operator-=(const BigInteger &a)
    {
        negative ^= 1;
        *this += a;
        negative ^= 1;
        normalize();
        return *this;
    }

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
};

const BigInteger operator+(const BigInteger &a, const BigInteger &b)
{
    BigInteger temp(a);
    temp += b;
    return temp;
}

const BigInteger operator-(const BigInteger &a, const BigInteger &b)
{
    BigInteger temp(a);
    temp -= b;
    return temp;
}

const BigInteger operator*(const BigInteger &a, const BigInteger &b)
{
    BigInteger temp(a);
    temp *= b;
    return temp;
}

const BigInteger operator/(const BigInteger &a, const BigInteger &b)
{
    BigInteger temp(a);
    temp /= b;
    return temp;
}

const BigInteger operator%(const BigInteger &a, const BigInteger &b)
{
    BigInteger temp(a);
    temp %= b;
    return temp;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &a)
{
    return out << a.toString();
}

std::istream &operator>>(std::istream &in, BigInteger &a)
{
    std::string s;
    in >> s;
    a = BigInteger(s);
    return in;
}

bool operator<(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) == CMP_LESS;
}

bool operator==(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) == CMP_EQUAL;
}

bool operator>(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) == CMP_GREATER;
}

bool operator<=(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) != CMP_GREATER;
}

bool operator>=(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) != CMP_LESS;
}

bool operator!=(const BigInteger &a, const BigInteger &b)
{
    return a._compare(b) != CMP_EQUAL;
}

const BigInteger abs(const BigInteger &a)
{
    return (a < 0) ? -a : a;
}

const BigInteger gcd(BigInteger a, BigInteger b)
{
    BigInteger d = 1;
    while (a && b)
    {
        if (!a.odd() && !b.odd())
        {
            d *= 2;
            a /= 2;
            b /= 2;
        }
        else if (!a.odd())
        {
            a /= 2;
        }
        else if (!b.odd())
        {
            b /= 2;
        }
        else if (a >= b)
        {
            a -= b;
        }
        else
        {
            b -= a;
        }
    }
    return d * (a ? a : b);
}

class Rational;
const Rational operator+(const Rational &a, const Rational &b);

class Rational
{
private:
    BigInteger _numerator, _denominator;

    void normalize()
    {
        if (_denominator < 0)
        {
            _denominator = -_denominator;
            _numerator = -_numerator;
        }
        BigInteger a = gcd(abs(_numerator), _denominator);
        _numerator /= a;
        _denominator /= a;
    }

    void expand(const BigInteger &a)
    {
        _numerator *= a;
        _denominator *= a;
    }

public:
    Rational(const BigInteger &a): _numerator(a), _denominator(1) {};

    Rational(long long a): _numerator(a), _denominator(1) {};

    Rational(): _numerator(0), _denominator(1) {};

    Rational(const BigInteger &a, const BigInteger &b): _numerator(a), _denominator(b)
    {
        normalize();
    }

    explicit Rational(const std::string &s): _numerator(), _denominator()
    {
        size_t i = 0;
        std::string s1, s2;
        for (; i < s.size() && (('0' <= s[i] && s[i] <= '9') || s[i] == '-'); ++i)
        {
            s1.push_back(s[i]);
        }
        if (i >= s.size())
        {
            *this = Rational(BigInteger(s1));
            return;
        }
        char c = s[i++];
        for (; i < s.size(); ++i)
        {
            s2.push_back(s[i]);
        }
        if (c == '/')
        {
            *this = Rational(BigInteger(s1), BigInteger(s2));
            return;
        }
        else
        {
            Rational fr = Rational(BigInteger(s2), BigInteger('1' + std::string(s2.size(), '0')));
            if (s1.size() && s1[0] == '-')
                fr = -fr;
            *this = Rational(BigInteger(s1)) + fr;
        }
    }

    compare_t compare(const Rational &a) const
    {
        return (_numerator * a._denominator - a._numerator * _denominator)._compare(0);
    }

    BigInteger &numerator()
    {
        return _numerator;
    }

    void numerator(const BigInteger &a)
    {
        _numerator = a;
        normalize();
    }

    BigInteger &denominator()
    {
        return _denominator;
    }

    void denominator(const BigInteger &a)
    {
        _denominator = a;
        normalize();
    }

    std::string toString() const
    {
        return _numerator.toString() + (_denominator == 1 ? "" : ("/" + _denominator.toString()));
    }

    explicit operator double() const
    {
        return double(_numerator) / double(_denominator);
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

    Rational &operator+=(const Rational &a)
    {
        BigInteger d = gcd(_denominator, a._denominator);
        expand(a._denominator / d);
        _numerator += a._numerator * _denominator / a._denominator;
        normalize();
        return *this;
    }

    Rational &operator-=(const Rational &a)
    {
        return *this += -a;
    }

    Rational &operator*=(const Rational &a)
    {
        _numerator *= a._numerator;
        _denominator *= a._denominator;
        normalize();
        return *this;
    }

    Rational &operator/=(const Rational &a)
    {
        if (!a._numerator)
            throw zero_division_error("Rational division by zero");
        _numerator *= a._denominator;
        _denominator *= a._numerator;
        normalize();
        return *this;
    }

    std::string asDecimal(size_t precision = 0) const
    {
        BigInteger newnum = _numerator;
        for (unsigned i = 0; i < precision; ++i)
        {
            newnum *= 10;
        }
        std::string result = abs(newnum / _denominator).toString();
        if (result.size() < precision + 1)
        {
            result = std::string(precision + 1 - result.size(), '0') + result;
        }
        if (precision)
            result.insert(result.end() - precision, '.');
        if (_numerator < 0)
            result = "-" + result;
        return result;
    }

};

bool operator<(const Rational &a, const Rational &b)
{
    return a.compare(b) == CMP_LESS;
}

bool operator==(const Rational &a, const Rational &b)
{
    return a.compare(b) == CMP_EQUAL;
}

bool operator>(const Rational &a, const Rational &b)
{
    return a.compare(b) == CMP_GREATER;
}

bool operator<=(const Rational &a, const Rational &b)
{
    return a.compare(b) != CMP_GREATER;
}

bool operator>=(const Rational &a, const Rational &b)
{
    return a.compare(b) != CMP_LESS;
}

bool operator!=(const Rational &a, const Rational &b)
{
    return a.compare(b) != CMP_EQUAL;
}

const Rational operator+(const Rational &a, const Rational &b)
{
    Rational temp(a);
    temp += b;
    return temp;
}

const Rational operator-(const Rational &a, const Rational &b)
{
    Rational temp(a);
    temp -= b;
    return temp;
}

const Rational operator*(const Rational &a, const Rational &b)
{
    Rational temp(a);
    temp *= b;
    return temp;
}

const Rational operator/(const Rational &a, const Rational &b)
{
    Rational temp(a);
    temp /= b;
    return temp;
}

std::ostream &operator<<(std::ostream &out, const Rational &a)
{
    return out << a.toString();
}

std::istream &operator>>(std::istream &in, Rational &a)
{
    std::string s;
    in >> s;
    a = Rational(s);
    return in;
}

#endif
