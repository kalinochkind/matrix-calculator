#include "biginteger.h"
#include <algorithm>

void BigInteger::normalize()
{
    while(digits.size() && digits.back() == 0)
        digits.pop_back();
    if(digits.empty())
        negative = false;
}

void BigInteger::shiftRight(size_t a)
{
    digits.resize(digits.size() + a);
    for(size_t i = digits.size() - 1; i >= a; --i)
    {
        digits[i] = digits[i - a];
    }
    for(size_t i = 0; i < a; i++)
    {
        digits[i] = 0;
    }
}

void BigInteger::multiply(unsigned a)
{
    long long carry = 0;
    for(size_t i = 0; carry || i < digits.size(); i++)
    {
        if(i >= digits.size())
            digits.push_back(0);
        long long new_carry = (digits[i] * 1ll * a + carry) / BLOCK_MOD;
        digits[i] = (digits[i] * 1ll * a + carry) % BLOCK_MOD;
        carry = new_carry;
    }
}

compare_t BigInteger::_compareAbs(const BigInteger &a) const
{
    if(digits.size() < a.digits.size())
        return CMP_LESS;
    if(digits.size() > a.digits.size())
        return CMP_GREATER;
    for(size_t i = digits.size(); i > 0; --i)
    {
        if(digits[i - 1] < a.digits[i - 1])
            return CMP_LESS;
        if(digits[i - 1] > a.digits[i - 1])
            return CMP_GREATER;
    }
    return CMP_EQUAL;
}

const BigInteger BigInteger::karatsuba(const BigInteger &a) const
{
    if(a.digits.size() <= 1)
    {
        BigInteger res(*this);
        res.multiply(a.digits[0]);
        return res;
    }
    size_t maxlen = std::max(digits.size(), a.digits.size());
    maxlen += maxlen % 2;
    BigInteger a1, a2, b1, b2;
    for(size_t i = 0; i < maxlen; i++)
    {
        if(i * 2 < maxlen)
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

void BigInteger::_divide(int a, int &remainder)
{
    if(!a)
        throw zero_division_error("Integer division by zero");
    if(digits.empty())
    {
        remainder = 0;
        return;
    }
    if(a < 0)
        negative ^= 1;
    long long carry = 0;
    for(auto i = digits.rbegin(); i != digits.rend(); ++i)
    {
        carry = carry * BLOCK_MOD + *i;
        *i = carry / abs(a);
        carry %= abs(a);
    }
    remainder = static_cast<int>(carry);
    normalize();
}

const BigInteger BigInteger::_divide_long(BigInteger a)
{
    negative = a.negative = false;
    BigInteger ans;
    if(_compareAbs(a) < 0)
        return 0;
    int shift = 0;
    while(_compareAbs(a) >= 0)
    {
        a.digits.insert(a.digits.begin(), 0);
        ++shift;
    }
    a.digits.erase(a.digits.begin());
    for(; shift; --shift)
    {
        int l = 0, r = BLOCK_MOD;
        while(l + 1 < r)
        {
            int m = (l + r) / 2;
            if(_compareAbs(a * m) < 0)
                r = m;
            else
                l = m;
        }
        ans.digits.push_back(l);
        *this -= a * l;
        a.digits.erase(a.digits.begin());
    }
    std::reverse(ans.digits.begin(), ans.digits.end());
    ans.normalize();
    return ans;
}

BigInteger::BigInteger(long long a): digits(), negative(false)
{
    if(a < 0)
    {
        negative = true;
        a = -a;
    }
    while(a)
    {
        digits.push_back(a % BLOCK_MOD);
        a /= BLOCK_MOD;
    }
}

BigInteger::BigInteger(const std::string &s): digits(), negative(false)
{
    unsigned temp = 0;
    unsigned pow = 1;
    for(auto i = s.crbegin(); i != s.crend(); ++i)
    {
        if(*i == '-' || *i == '+')
        {
            if(*i == '-')
                negative = true;
            continue;
        }
        if(*i < '0' || *i > '9')
            throw invalid_number_error("Invalid integer: " + s);
        temp += pow * (*i - '0');
        pow *= 10;
        if(pow == BLOCK_MOD)
        {
            pow = 1;
            digits.push_back(temp);
            temp = 0;
        }
    }
    if(s == "-")
        temp = 1;
    digits.push_back(temp);
    normalize();
}

compare_t BigInteger::_compare(const BigInteger &a) const
{
    if(negative && !a.negative)
        return CMP_LESS;
    if(!negative && a.negative)
        return CMP_GREATER;
    return static_cast<compare_t>(_compareAbs(a) * (negative ? -1 : 1));
}

std::string BigInteger::toString() const
{
    if(digits.empty())
        return "0";
    std::string res = (negative ? "-" : "");
    for(auto i = digits.rbegin(); i != digits.rend(); ++i)
    {
        std::string temp = ::toString(*i);
        if(i != digits.rbegin())
        {
            while(temp.length() < BLOCK_SIZE)
            {
                temp = "0" + temp;
            }
        }
        res += temp;
    }
    return res;
}

BigInteger &BigInteger::operator+=(const BigInteger &a)
{
    bool sub = negative ^a.negative;
    int rem = 0;
    if(!sub)  // addition
    {
        for(size_t i = 0; i < std::max(a.digits.size(), digits.size()) || rem; ++i)
        {
            if(digits.size() <= i)
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
        for(size_t i = 0; i < a.digits.size() || rem; ++i)
        {
            if(digits.size() <= i)
                digits.push_back(0);
            unsigned dig = (small.digits.size() > i ? small.digits[i] : 0);
            if(big.digits[i] >= dig + rem)
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

BigInteger &BigInteger::operator*=(BigInteger a)
{
    bool neg = negative ^a.negative;
    if(_compareAbs(a) == CMP_LESS)
    {
        std::swap(*this, a);
    }
    if(!a.digits.size())
    {
        digits.clear();
    }
    else if(a.digits.size() == 1)
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

BigInteger &BigInteger::operator/=(const BigInteger &a)
{
    if(!a)
        throw zero_division_error("Integer division by zero");
    if(a.digits.size() <= 1)
    {
        int temp;
        _divide(a.digits[0], temp);
        return *this;
    }
    bool neg = a.negative ^negative;
    *this = _divide_long(a);
    negative = neg;
    normalize();
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &a)
{
    BigInteger temp(*this);
    temp /= a;
    temp *= a;
    return *this -= temp;
}

BigInteger &BigInteger::operator-=(const BigInteger &a)
{
    negative ^= 1;
    *this += a;
    negative ^= 1;
    normalize();
    return *this;
}

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


const BigInteger gcd(BigInteger a, BigInteger b)
{
    BigInteger d = 1;
    while(a && b)
    {
        if(!a.odd() && !b.odd())
        {
            d *= 2;
            a /= 2;
            b /= 2;
        }
        else if(!a.odd())
        {
            a /= 2;
        }
        else if(!b.odd())
        {
            b /= 2;
        }
        else if(a >= b)
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


std::pair<BigInteger, BigInteger> ext_gcd(BigInteger a, BigInteger b)
{
    BigInteger u = 1, v = 0, s = 0, t = 1;
    while(!a.odd() && !b.odd())
    {
        a /= 2;
        b /= 2;
    }
    BigInteger alpha = a, beta = b;
    while(!a.odd())
    {
        a /= 2;
        if(!u.odd() && !b.odd())
        {
            u /= 2;
            v /= 2;
        }
        else
        {
            u += beta;
            u /= 2;
            v -= alpha;
            v /= 2;
        }
    }
    while(a != b)
    {
        if(!b.odd())
        {
            b /= 2;
            if(!s.odd() && !t.odd())
            {
                s /= 2;
                t /= 2;
            }
            else
            {
                s += beta;
                s /= 2;
                t -= alpha;
                t /= 2;
            }
        }
        else if(b < a)
        {
            a.swap(b);
            u.swap(s);
            v.swap(t);
        }
        else
        {
            b -= a;
            s -= u;
            t -= v;
        }
    }
    return {s, t};
}