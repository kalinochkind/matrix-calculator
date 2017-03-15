#include "rational.h"

void Rational::normalize()
{
    if(_denominator < 0)
    {
        _denominator = -_denominator;
        _numerator = -_numerator;
    }
    BigInteger a = gcd(abs(_numerator), _denominator);
    _numerator /= a;
    _denominator /= a;
}

void Rational::expand(const BigInteger &a)
{
    _numerator *= a;
    _denominator *= a;
}

Rational::Rational(const std::string &s): _numerator(), _denominator()
{
    size_t i = 0;
    std::string s1, s2;
    for(; i < s.size() && (('0' <= s[i] && s[i] <= '9') || s[i] == '-' || s[i] == '+'); ++i)
    {
        s1.push_back(s[i]);
    }
    if(i >= s.size())
    {
        *this = Rational(BigInteger(s1));
        return;
    }
    char c = s[i++];
    for(; i < s.size(); ++i)
    {
        if(s[i] < '0' || s[i] > '9')
            throw invalid_number_error("Invalid rational: " + s);
        s2.push_back(s[i]);
    }
    if(c == '/')
    {
        *this = Rational(BigInteger(s1), BigInteger(s2));
        return;
    }
    else if(c == '.')
    {
        Rational fr = Rational(BigInteger(s2), BigInteger('1' + std::string(s2.size(), '0')));
        if(s1.size() && s1[0] == '-')
            fr = -fr;
        *this = Rational(BigInteger(s1)) + fr;
    }
    else
        throw invalid_number_error("Invalid rational: " + s);
}

std::string Rational::toString() const
{
    return _numerator.toString() + (_denominator == 1 ? "" : ("/" + _denominator.toString()));
}

Rational &Rational::operator+=(const Rational &a)
{
    BigInteger d = gcd(_denominator, a._denominator);
    expand(a._denominator / d);
    _numerator += a._numerator * _denominator / a._denominator;
    normalize();
    return *this;
}

Rational &Rational::operator*=(const Rational &a)
{
    _numerator *= a._numerator;
    _denominator *= a._denominator;
    if(_denominator != 1)
        normalize();
    return *this;
}

Rational &Rational::operator/=(const Rational &a)
{
    if(!a._numerator)
        throw zero_division_error("Rational division by zero");
    _numerator *= a._denominator;
    _denominator *= a._numerator;
    normalize();
    return *this;
}

const std::string Rational::asDecimal(size_t precision) const
{
    BigInteger newnum = _numerator;
    for(unsigned i = 0; i < precision; ++i)
    {
        newnum *= 10;
    }
    std::string result = abs(newnum / _denominator).toString();
    if(result.size() < precision + 1)
    {
        result = std::string(precision + 1 - result.size(), '0') + result;
    }
    if(precision)
        result.insert(result.end() - precision, '.');
    if(result.size() > 1 && result != "0." + std::string(result.size() - 2, '0'))
        while(precision && result.back() == '0')
            result.pop_back();
    if(_numerator < 0)
        result = "-" + result;
    return result;
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