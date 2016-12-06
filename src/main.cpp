#include "matrix.h"
#include "parser.h"
#include "finite.h"
#include <iostream>
#include <map>
#include <sstream>
#include <cstdlib>

using namespace std;

void die(const string &s)
{
    cout << s << endl;
    exit(1);
}

string safeGetline()
{
    string s;
    while (s.empty())
    {
        getline(cin, s);
        while (s.size() && isspace(s.back()))
            s.pop_back();
    }
    return s;
}

template<class FMatrix>
FMatrix getMatrix(string prompt)
{
    string s, sum;
    cout << prompt << endl;
    s = safeGetline();
    unsigned width = 0, height = 0;
    while (s.length())
    {
        istringstream is;
        is.str(s);
        unsigned cwidth = 0;
        Rational dummy;
        while (is >> dummy)
            ++cwidth;
        if (!cwidth)
            break;
        if (width && width != cwidth)
            die("Incorrect matrix");
        width = cwidth;
        ++height;
        sum += s + ' ';
        getline(cin, s);
    }
    FMatrix m(height, width);
    istringstream is;
    is.str(sum);
    is >> m;
    return m;
}

template<class NumMatrix>
void processOp(string op, vector<NumMatrix> &st,
               map<string, pair<int, NumMatrix (*)(const vector<NumMatrix *> &)> > &operations)
{
    if (operations[op].first == 1)
    {
        st.back() = operations[op].second({&st.back()});
    }
    else if (operations[op].first == 2)
    {
        NumMatrix a = st.back();
        st.pop_back();
        st.back() = operations[op].second({&st.back(), &a});
    }
    else
    {
        NumMatrix a = st.back();
        st.pop_back();
        NumMatrix b = st.back();
        st.pop_back();
        st.back() = operations[op].second({&st.back(), &b, &a});
    }
}

template<class Field>
struct _NumMatrix
{
    int im;
    Matrix<Field> fm;
    bool is_int;

    _NumMatrix(): im(0), fm(), is_int(true) {}

    _NumMatrix(const _NumMatrix &a): im(a.im), fm(a.fm), is_int(a.is_int) {}

    _NumMatrix(int a): im(a), fm(), is_int(true) {}

    _NumMatrix(const Matrix<Field> &a): im(0), fm(a), is_int(false) {}

    _NumMatrix(const Field &a): im(0), fm(Matrix<Field>::fromNumber(a)), is_int(false) {}

    _NumMatrix &operator=(const _NumMatrix &a)
    {
        is_int = a.is_int;
        im = a.im;
        fm = a.fm;
        return *this;
    }

    Matrix<Field> toMatrix() const
    {
        if (is_int)
            return Matrix<Field>::fromNumber(Field(im));
        else
            return fm;
    }

    _NumMatrix operator*(const _NumMatrix &m) const
    {
        if (is_int && m.is_int)
            return im * m.im;
        if (is_int)
            return m.fm * Field(im);
        if (m.is_int)
            return fm * Field(m.im);
        if (fm.width() == 1 && fm.height() == 1)
            return m.fm * fm[0][0];
        if (m.fm.width() == 1 && m.fm.height() == 1)
            return fm * m.fm[0][0];
        return fm * m.fm;
    }

    _NumMatrix operator/(const _NumMatrix &m) const
    {
        return operator*(m.toMatrix().inverted());
    }

    _NumMatrix operator+(const _NumMatrix &m) const
    {
        if (is_int && m.is_int)
            return im + m.im;
        return toMatrix() + m.toMatrix();
    }

    _NumMatrix operator-() const
    {
        if (is_int)
            return -im;
        else
            return -fm;
    }

};

Matrix<Finite> f_cfrac(Finite a)
{
    return Matrix<Finite>::fromNumber(a);
}

Matrix<Rational> f_cfrac(Rational a)
{
    std::vector<Rational> ans = {0};
    if (a < 0)
    {
        ans[0] -= (-a.numerator() / a.denominator()) + 1;
        a += (-a.numerator() / a.denominator()) + 1;
    }
    ans[0] += (a.numerator() / a.denominator());
    a -= (a.numerator() / a.denominator());
    while (a.denominator() != 1)
    {
        a = 1 / a;
        ans.push_back(a.numerator() / a.denominator());
        a -= (a.numerator() / a.denominator());
    }
    Matrix<Rational> res(1, ans.size());
    for (unsigned i = 0; i < ans.size(); ++i)
    {
        res[0][i] = ans[i];
    }
    return res;
}


template<class T>
T f_revcfrac(const Matrix<T> &a)
{
    T res = a[0][a.width() - 1];
    for (int i = int(a.width()) - 2; i >= 0; --i)
    {
        res = T(1) / res + a[0][i];
    }
    return res;
}


template<class Field>
void f_expr()
{
    typedef _NumMatrix<Field> NumMatrix;
    map<string, pair<int, NumMatrix (*)(const vector<NumMatrix *> &)> > operations =
            {
                    {"+",         {2, [](const vector<NumMatrix *> &a) { return *a[0] + *a[1]; }}},
                    {"^",         {2, [](const vector<NumMatrix *> &a) {
                        if (!a[1]->is_int)
                            die("Invalid use of ^: integer required");
                        return NumMatrix(a[0]->toMatrix().power(a[1]->im));
                    }}},
                    {"*",         {2, [](const vector<NumMatrix *> &a) {
                        return *a[0] * *a[1];
                    }}},
                    {"/",         {2, [](const vector<NumMatrix *> &a) {
                        return *a[0] / *a[1];
                    }}},
                    {"-",         {2, [](const vector<NumMatrix *> &a) { return *a[0] + -*a[1]; }}},
                    {"_",         {1, [](const vector<NumMatrix *> &a) { return -*a[0]; }}},
                    {"det",       {1, [](const vector<NumMatrix *> &a) { return NumMatrix(a[0]->toMatrix().det()); }}},
                    {"rank",      {1, [](const vector<NumMatrix *> &a) { return NumMatrix(a[0]->toMatrix().rank()); }}},
                    {"trace",     {1, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().trace());
                    }}},
                    {"transpose", {1, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().transposed());
                    }}},
                    {"id",        {1, [](const vector<NumMatrix *> &a) {
                        if (!a[1]->is_int)
                            die("Invalid use of id");
                        return NumMatrix(Matrix<Field>::identity(abs(a[0]->im)));
                    }}},
                    {"=",         {2, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix() == a[1]->toMatrix());
                    }}},
                    {"width",     {1, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().width());
                    }}},
                    {"height",    {1, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().height());
                    }}},
                    {"solve",     {1, [](const vector<NumMatrix *> &a) {
                        unsigned sz = a[0]->toMatrix().height();
                        if (!sz || a[0]->toMatrix().width() != sz + 1)
                            die("Invalid use of solve: N*N+1 matrix required");
                        Matrix<Field> sys = a[0]->toMatrix().submatrix(0, 0, sz - 1, sz - 1);
                        Matrix<Field> right = a[0]->toMatrix().submatrix(0, sz, sz - 1, sz);
                        sys.inverseExt(right);
                        return NumMatrix(right.transposed());
                    }}},
                    {"at",        {3, [](const vector<NumMatrix *> &a) {
                        if (!a[1]->is_int || !a[2]->is_int)
                            die("Invalid use of at");
                        if (a[1]->im < 0 || a[1]->im >= int(a[0]->toMatrix().height()) ||
                            a[2]->im < 0 || a[2]->im >= int(a[0]->toMatrix().width()))
                            die("at: out of range");
                        return NumMatrix(a[0]->toMatrix()[a[1]->im][a[2]->im]);
                    }}},
                    {"int",       {1, [](const vector<NumMatrix *> &a) {
                        if (a[0]->is_int)
                            return *a[0];
                        if (a[0]->fm.width() != 1 || a[0]->fm.height() != 1)
                            die("int: matrix 1*1 required");
                        return NumMatrix(int(a[0]->fm[0][0]));
                    }}},
                    {"cfrac",     {1, [](const vector<NumMatrix *> &a) {
                        auto m = a[0]->toMatrix();
                        if (m.width() != 1 || m.height() != 1)
                            die("cfrac: matrix 1*1 required");
                        return NumMatrix(f_cfrac(m[0][0]));
                    }}},
                    {"rcfrac",    {1, [](const vector<NumMatrix *> &a) {
                        auto m = a[0]->toMatrix();
                        if (m.height() != 1 || !m.width())
                            die("revcfrac: matrix 1*n required");
                        return NumMatrix(f_revcfrac(m));
                    }}},
                    {"joinh",     {2, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().joinHorizontal(a[1]->toMatrix()));
                    }}},
                    {"joinv",     {2, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().joinVertical(a[1]->toMatrix()));
                    }}},
            };
    cout << "Expression: ";
    string s = safeGetline();
    auto v = splitExpression(s);
    map<char, NumMatrix> mmap;
    vector<pair<token_type, string> > opst;
    vector<NumMatrix> st;
    int st_size = 0;
    for (pair<token_type, string> &i : v)
    {
        switch (i.first)
        {
            case TOKEN_NUMBER:
            case TOKEN_MATRIX:
                ++st_size;
                break;
            case TOKEN_OP:
                while (opst.size() && opst.back().first == TOKEN_OP &&
                       priority[int(i.second[0])] + rightassoc[int(i.second[0])] <=
                       priority[int(opst.back().second[0])])
                {
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
                }
                if (st_size < 0)
                    die("Invalid expression");
            case TOKEN_FUNC:
                if (!operations.count(i.second))
                    die("Invalid function: " + i.second);
            case TOKEN_LEFTPAR:
                opst.push_back(i);
                break;
            case TOKEN_RIGHTPAR:
                while (opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
                }
                if (opst.empty() || st_size <= 0)
                    die("Invalid expression");
                opst.pop_back();
                if (opst.size() && opst.back().first == TOKEN_FUNC)
                {
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
                }
                break;
            case TOKEN_COMMA:
                while (opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
                }
                if (opst.empty() || st_size <= 0)
                    die("Invalid expression");
                break;
        }
    }
    while (opst.size())
    {
        if (opst.back().first == TOKEN_LEFTPAR || opst.back().first == TOKEN_RIGHTPAR)
            die("Invalid expression");
        st_size -= operations[opst.back().second].first - 1;
        opst.pop_back();
    }
    if (st_size != 1)
        die("Invalid expression");
    for (pair<token_type, string> &i : v)
    {
        Rational tt;
        istringstream is, iis;
        switch (i.first)
        {
            case TOKEN_NUMBER:
                is.str(i.second);
                is >> tt;
                if (tt == int(tt))
                    st.push_back(NumMatrix(int(tt)));
                else
                {
                    Field tf;
                    iis.str(i.second);
                    iis >> tf;
                    st.push_back(NumMatrix(tf));
                }
                break;
            case TOKEN_MATRIX:
                if (!mmap.count(i.second[0]))
                    mmap[i.second[0]] = getMatrix<Matrix<Field>>(string("Matrix ") + i.second + ':');
                st.push_back(mmap[i.second[0]]);
                break;
            case TOKEN_OP:
                while (opst.size() && opst.back().first == TOKEN_OP &&
                       priority[int(i.second[0])] + rightassoc[int(i.second[0])] <=
                       priority[int(opst.back().second[0])])
                {
                    processOp(opst.back().second, st, operations);
                    opst.pop_back();
                }
            case TOKEN_FUNC:
            case TOKEN_LEFTPAR:
                opst.push_back(i);
                break;
            case TOKEN_RIGHTPAR:
                while (opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    processOp(opst.back().second, st, operations);
                    opst.pop_back();
                }
                opst.pop_back();
                if (opst.size() && opst.back().first == TOKEN_FUNC)
                {
                    processOp(opst.back().second, st, operations);
                    opst.pop_back();
                }
                break;
            case TOKEN_COMMA:
                while (opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    processOp(opst.back().second, st, operations);
                    opst.pop_back();
                }
                break;
        }
    }
    while (opst.size())
    {
        processOp(opst.back().second, st, operations);
        opst.pop_back();
    }
    cout << "Result:\n" << st[0].toMatrix();
}


int main(int argc, char **argv)
{
    try
    {
        if (argc == 1)
        {
            f_expr<Rational>();
        }
        else
        {
            _FINITE_ORDER = atoi(argv[1]);
            if (_FINITE_ORDER < 2)
                die("Order must be at least 2");
            f_expr<Finite>();
        }
    }
    catch (matrix_error e)
    {
        cout << "Matrix error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
