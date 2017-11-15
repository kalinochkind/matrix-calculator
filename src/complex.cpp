#include "complex.h"

const Complex Complex::inverse() const
{
    Rational sqmod = _re * _re + _im * _im;
    if(!sqmod)
        throw zero_division_error("Complex division by zero");
    return Complex(_re / sqmod, -_im / sqmod);
}


std::string Complex::toString() const
{
    if(!_im)
    {
        return _re.toString();
    }
    if(!_re)
    {
        if(_im == 1)
            return "I";
        if(_im == -1)
            return "-I";
        return _im.toString() + "I";
    }
    std::string s = _re.toString() + (_im >= 0 ? "+" : "");
    if(_im == 1)
        return s + "I";
    if(_im == -1)
        return s + "-I";
    return s + _im.toString() + "I";
}


Complex::Complex(const std::string &s): _re(), _im()
{
    size_t i = 0;
    std::string s1, s2;
    for(; i < s.size() && (i == 0 || (s[i] != '+' && s[i] != '-')) && s[i] != 'I'; ++i)
    {
        s1.push_back(s[i]);
    }
    if(i >= s.size())
    {
        *this = Complex(Rational(s1));
        return;
    }
    if(s[i] == 'I' && i + 1 == s.size())
    {
        if(s1.empty())
            s1 = "1";
        *this = Complex(0, Rational(s1));
        return;
    }
    for(; i < s.size() && s[i] != 'I'; ++i)
    {
        s2.push_back(s[i]);
    }
    if(s[i] != 'I' || i + 1 != s.size())
    {
        throw invalid_number_error("Invalid complex: " + s);
    }
    *this = Complex(Rational(s1), Rational(s2));
}

Complex &Complex::operator+=(const Complex &a)
{
    _re += a._re;
    _im += a._im;
    return *this;
}

Complex &Complex::operator-=(const Complex &a)
{
    _re -= a._re;
    _im -= a._im;
    return *this;
}

Complex &Complex::operator*=(const Complex &a)
{
    Rational new_re = _re * a._re - _im * a._im;
    _im = _re * a._im + _im * a._re;
    _re = std::move(new_re);
    return *this;
}

Complex &Complex::operator/=(const Complex &a)
{
    return *this *= a.inverse();
}

const Complex operator+(const Complex &a, const Complex &b)
{
    Complex t = a;
    return t += b;
}

const Complex operator-(const Complex &a, const Complex &b)
{
    Complex t = a;
    return t -= b;
}

const Complex operator*(const Complex &a, const Complex &b)
{
    Complex t = a;
    return t *= b;
}

const Complex operator/(const Complex &a, const Complex &b)
{
    Complex t = a;
    return t /= b;
}

std::ostream &operator<<(std::ostream &out, const Complex &a)
{
    return out << a.toString();
}

std::istream &operator>>(std::istream &in, Complex &a)
{
    std::string s;
    in >> s;
    a = Complex(s);
    return in;
}