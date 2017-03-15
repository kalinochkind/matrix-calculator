#include "matrix.h"
#include "parser.h"
#include "finite.h"
#include "polynom.h"
#include <map>
#include <set>
#include <sstream>
#include <cstdlib>
#include <readline/readline.h>
#include <readline/history.h>

using namespace std;

BigInteger _FINITE_ORDER = 2;

class end_of_input: public exception
{
    using exception::exception;
};

class invalid_expression: public exception
{
    using exception::exception;
};

void die(const string &s)
{
    cout << s << endl;
    throw invalid_expression();
}

bool isdigit(const string &s)
{
    for(char i : s)
    {
        if(i < '0' || i > '9')
            return false;
    }
    return not s.empty();
}

string _readline(const string &prompt = "")
{
    char *cs = readline(prompt.c_str());
    if(!cs)
    {
        cout << endl;
        throw end_of_input();
    }
    string s = cs;
    free(cs);
    return s;
}

string safeGetline(const string &prompt = "")
{
    string s;
    while(s.empty())
    {
        s = _readline(prompt);
        while(s.size() && isspace(s.back()))
            s.pop_back();
    }
    return s;
}


template<class Field>
Matrix<Field> getMatrix(string prompt)
{
    string s, sum;
    cout << prompt << endl;
    s = safeGetline();
    unsigned width = 0, height = 0;
    while(s.length())
    {
        istringstream is;
        is.str(s);
        unsigned cwidth = 0;
        Field dummy;
        while(is >> dummy)
            ++cwidth;
        if(!cwidth)
            break;
        if(width && width != cwidth)
            throw matrix_error("Incorrect matrix");
        width = cwidth;
        ++height;
        sum += s + ' ';
        s = _readline();
    }
    Matrix<Field> m(height, width);
    istringstream is;
    is.str(sum);
    is >> m;
    return m;
}

template<class NumMatrix>
void processOp(string op, vector<NumMatrix> &st, size_t opcount,
               map<string, pair<int, NumMatrix (*)(const vector<NumMatrix *> &)> > &operations)
{
    if(!opcount)
    {
        if(operations[op].first <= 0)
            die("Invalid number of arguments: " + op);
        opcount = operations[op].first;
    }
    vector<NumMatrix*> args;
    for(;opcount;--opcount)
    {
        args.push_back(&st[st.size() - opcount]);
    }
    auto res = operations[op].second(args);
    for(size_t i=0;i<args.size();++i)
    {
        st.pop_back();
    }
    st.push_back(res);
}

enum class NumMatrixType
{
    number, matrix, polynom
};

template<class Field>
struct _NumMatrix
{
    BigInteger im;
    Matrix<Field> fm;
    Polynom<Field> pm;
    NumMatrixType type;

    _NumMatrix(): im(0), fm(), pm(), type(NumMatrixType::number) {}

    _NumMatrix(const _NumMatrix &a): im(a.im), fm(a.fm), pm(a.pm), type(a.type) {}

    _NumMatrix(const BigInteger &a): im(a), fm(), pm(), type(NumMatrixType::number) {}

    _NumMatrix(unsigned a): im(a), fm(), pm(), type(NumMatrixType::number) {}

    _NumMatrix(int a): im(a), fm(), pm(), type(NumMatrixType::number) {}

    _NumMatrix(const Matrix<Field> &a): im(0), fm(a), pm(), type(NumMatrixType::matrix) {}

    _NumMatrix(const Polynom<Field> &a): im(0), fm(), pm(a), type(NumMatrixType::polynom) {}

    _NumMatrix(const Field &a): im(0), fm(Matrix<Field>::fromNumber(a)), pm(), type(NumMatrixType::matrix) {}

    _NumMatrix &operator=(const _NumMatrix &a)
    {
        type = a.type;
        im = a.im;
        fm = a.fm;
        pm = a.pm;
        return *this;
    }

    Matrix<Field> toMatrix() const
    {
        if(type == NumMatrixType::number)
            return Matrix<Field>::fromNumber(Field(im));
        else if(type == NumMatrixType::polynom)
            return pm.toMatrix();
        else
            return fm;
    }

    const _NumMatrix operator*(const _NumMatrix &m) const
    {
        switch(type)
        {
            case NumMatrixType::number:
                if(m.type == NumMatrixType::number)
                    return im * m.im;
                if(m.type == NumMatrixType::matrix)
                    return m.fm * Field(im);
                return m.pm * Field(im);
            case NumMatrixType::matrix:
                if(m.type == NumMatrixType::number)
                    return fm * Field(m.im);
                if(m.type == NumMatrixType::matrix)
                {
                    if(fm.width() == 1 && fm.height() == 1)
                        return m.fm * fm[0][0];
                    else if(m.fm.width() == 1 && m.fm.height() == 1)
                        return fm * m.fm[0][0];
                    else
                        return fm * m.fm;
                }
                if(fm.width() == 1 && fm.height() == 1)
                    return m.pm * fm[0][0];
                if(pm.degree() < 1)
                    return fm * m.pm.toMatrix()[0][0];
                die("Matrix * polynom is undefined");
            case NumMatrixType::polynom:
                if(m.type == NumMatrixType::number)
                    return pm * Field(m.im);
                if(m.type == NumMatrixType::matrix)
                {
                    if(m.fm.width() == 1 && m.fm.height() == 1)
                        return pm * m.fm[0][0];
                    else if(pm.degree() < 1)
                        return m.fm * pm.toMatrix()[0][0];
                    else
                        die("Polynom * matrix is undefined");
                }
                return pm * m.pm;
        }
        return 0;
    }

    const _NumMatrix operator/(const _NumMatrix &m) const
    {
        if(m.type == NumMatrixType::number || m.type == NumMatrixType::matrix)
            return operator*(m.toMatrix().inverted());
        if(type == NumMatrixType::polynom)
            return pm / m.pm;
        if(type == NumMatrixType::matrix && (fm.width() != 1 || fm.height() != 1))
            die("Polynom / matrix is undefined");
        return pm * m.toMatrix().inverted()[0][0];
    }

    const _NumMatrix operator%(const _NumMatrix &m) const
    {
        if(m.type == NumMatrixType::number || m.type == NumMatrixType::matrix)
            return im % m.im;
        if(type == NumMatrixType::polynom && m.type == NumMatrixType::polynom)
            return pm % m.pm;
        die("Invalid use of %");
        return 0;
    }

    const _NumMatrix operator+(const _NumMatrix &m) const
    {
        if(type == NumMatrixType::number && m.type == NumMatrixType::number)
            return im + m.im;
        if(type == NumMatrixType::polynom && m.type == NumMatrixType::polynom)
            return pm + m.pm;
        Matrix<Field> m1 = toMatrix(), m2 = m.toMatrix();
        if(m1.width() == 1 && m1.height() == 1 && m2.width() == m2.height())
            return m2 + Matrix<Field>::identity(m2.width()) * m1[0][0];
        if(m2.width() == 1 && m2.height() == 1 && m1.width() == m1.height())
            return m1 + Matrix<Field>::identity(m1.width()) * m2[0][0];
        return m1 + m2;
    }

    const _NumMatrix operator-() const
    {
        if(type == NumMatrixType::number)
            return -im;
        else if(type == NumMatrixType::matrix)
            return -fm;
        else
            return -pm;
    }

};

Matrix<Finite> f_cfrac(Finite)
{
    die("cfrac is supported only for rationals");
    return 0;
}

Matrix<Rational> f_cfrac(Rational a)
{
    std::vector<Rational> ans = {0};
    if(a < 0)
    {
        ans[0] -= (-a.numerator() / a.denominator()) + 1;
        a += (-a.numerator() / a.denominator()) + 1;
    }
    ans[0] += (a.numerator() / a.denominator());
    a -= (a.numerator() / a.denominator());
    while(a.denominator() != 1)
    {
        a = 1 / a;
        ans.push_back(a.numerator() / a.denominator());
        a -= (a.numerator() / a.denominator());
    }
    Matrix<Rational> res(1, ans.size());
    for(unsigned i = 0; i < ans.size(); ++i)
    {
        res[0][i] = ans[i];
    }
    return res;
}


Rational f_revcfrac(const Matrix<Rational> &a)
{
    Rational res = a[0][a.width() - 1];
    for(int i = int(a.width()) - 2; i >= 0; --i)
    {
        res = 1 / res + a[0][i];
    }
    return res;
}

Finite f_revcfrac(const Matrix<Finite> &)
{
    die("rcfrac is supported only for rationals");
    return 0;
}

void printDecimalResult(const Finite &) {}

void printDecimalResult(const Rational &a)
{
    if(a.denominator() == 1)
        return;
    cout << "decimal: " << a.asDecimal(30) << endl;
}

template<class Field>
Polynom<Field> gcd2(const Polynom<Field> &a, const Polynom<Field> &b)
{
    if(a.degree() < 1 && b.degree() < 1)
    {
        BigInteger i1(a.toMatrix()[0][0]);
        BigInteger i2(b.toMatrix()[0][0]);
        return Polynom<Field>(Field(gcd(i1, i2)));
    }
    return a.gcd(b);
}

template<class Field>
void f_expr(string expr)
{
    typedef _NumMatrix<Field> NumMatrix;
    map<string, pair<int, NumMatrix (*)(const vector<NumMatrix *> &)> > operations =
            {
                    {"+",         {2, [](const vector<NumMatrix *> &a) { return *a[0] + *a[1]; }}},
                    {"^",         {2, [](const vector<NumMatrix *> &a) {
                        if(a[0]->type == NumMatrixType::polynom || a[1]->type == NumMatrixType::polynom)
                            die("^ is not supported for polynoms");
                        if(a[1]->type == NumMatrixType::number)
                            return NumMatrix(a[0]->toMatrix().power(a[1]->im));
                        Matrix<Field> m = a[1]->toMatrix();
                        if(m.height() != 1 || m.width() != 1)
                            die("Invalid use of ^: integer required");
                        return NumMatrix(a[0]->toMatrix().power(BigInteger(m[0][0])));
                    }}},
                    {"*",         {2, [](const vector<NumMatrix *> &a) {
                        return *a[0] * *a[1];
                    }}},
                    {"/",         {2, [](const vector<NumMatrix *> &a) {
                        return *a[0] / *a[1];
                    }}},
                    {"%",         {2, [](const vector<NumMatrix *> &a) {
                        return *a[0] % *a[1];
                    }}},
                    {"-",         {2, [](const vector<NumMatrix *> &a) { return *a[0] + -*a[1]; }}},
                    {"\\",        {1, [](const vector<NumMatrix *> &a) { return -*a[0]; }}},
                    {"det",       {1, [](const vector<NumMatrix *> &a) { return NumMatrix(a[0]->toMatrix().det()); }}},
                    {"rank",      {1, [](const vector<NumMatrix *> &a) { return NumMatrix(a[0]->toMatrix().rank()); }}},
                    {"trace",     {1, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().trace());
                    }}},
                    {"diag",      {-1, [](const vector<NumMatrix *> &a) {
                        Matrix<Field> m(a.size());
                        for(unsigned i=0;i<a.size();++i)
                        {
                            m[i][i] = a[i]->toMatrix()[0][0];
                        }
                        return NumMatrix(m);
                    }}},
                    {"transpose", {1, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().transposed());
                    }}},
                    {"id",        {1, [](const vector<NumMatrix *> &a) {
                        if(a[0]->type != NumMatrixType::number || a[0]->im <= 0)
                            die("Invalid use of id");
                        return NumMatrix(Matrix<Field>::identity(int(a[0]->im)));
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
                        return NumMatrix(a[0]->toMatrix().partial());
                    }}},
                    {"at",        {3, [](const vector<NumMatrix *> &a) {
                        if(a[1]->type != NumMatrixType::number || a[2]->type != NumMatrixType::number)
                            die("Invalid use of at");
                        if(a[1]->im < 0 || a[1]->im >= int(a[0]->toMatrix().height()) ||
                           a[2]->im < 0 || a[2]->im >= int(a[0]->toMatrix().width()))
                            die("at: out of range");
                        return NumMatrix(a[0]->toMatrix()[int(a[1]->im)][int(a[2]->im)]);
                    }}},
                    {"int",       {1, [](const vector<NumMatrix *> &a) {
                        if(a[0]->type == NumMatrixType::number)
                            return *a[0];
                        if(a[0]->type == NumMatrixType::polynom)
                            die("at: number or matrix 1*1 required");
                        if(a[0]->fm.width() != 1 || a[0]->fm.height() != 1)
                            die("int: matrix 1*1 required");
                        return NumMatrix(int(a[0]->fm[0][0]));
                    }}},
                    {"cfrac",     {1, [](const vector<NumMatrix *> &a) {
                        auto m = a[0]->toMatrix();
                        if(m.width() != 1 || m.height() != 1)
                            die("cfrac: matrix 1*1 required");
                        return NumMatrix(f_cfrac(m[0][0]));
                    }}},
                    {"rcfrac",    {1, [](const vector<NumMatrix *> &a) {
                        auto m = a[0]->toMatrix();
                        if(m.height() != 1 || !m.width())
                            die("rcfrac: matrix 1*n required");
                        return NumMatrix(f_revcfrac(m));
                    }}},
                    {"joinh",     {-1, [](const vector<NumMatrix *> &a) {
                        Matrix<Field> res = a[0]->toMatrix();
                        for(size_t i=1;i<a.size();++i)
                            res = res.joinHorizontal(a[i]->toMatrix());
                        return NumMatrix(res);
                    }}},
                    {"joinv",     {-1, [](const vector<NumMatrix *> &a) {
                        Matrix<Field> res = a[0]->toMatrix();
                        for(size_t i=1;i<a.size();++i)
                            res = res.joinVertical(a[i]->toMatrix());
                        return NumMatrix(res);
                    }}},
                    {"gauss",     {1, [](const vector<NumMatrix *> &a) {
                        Matrix<Field> m = a[0]->toMatrix();
                        m.gauss();
                        return NumMatrix(m);
                    }}},
                    {"fund",     {1, [](const vector<NumMatrix *> &a) {
                        Matrix<Field> m = a[0]->toMatrix();
                        return NumMatrix(m.fundamental());
                    }}},
                    {"gcd",       {-1, [](const vector<NumMatrix *> &a) {
                        if(a.size() < 2)
                            die("gcd requires at least 2 arguments");
                        Polynom<Field> p(a[0]->toMatrix());
                        for(size_t i=1;i<a.size();++i)
                            p = gcd2(p, Polynom<Field>(a[i]->toMatrix()));
                        return NumMatrix(p);
                    }}},
                    {"degree",    {1, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix());
                        return NumMatrix(p1.degree());
                    }}},
                    {"polynom",   {1, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix());
                        return NumMatrix(p1);
                    }}},
                    {"divmod", {2, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix()), p2(a[1]->toMatrix());
                        Polynom<Field> div = p1 / p2, mod = p1 % p2;
                        unsigned mx = std::max(div.degree(), mod.degree()) + 1;
                        return NumMatrix(div.toMatrix(mx).joinVertical(mod.toMatrix(mx)));
                    }}},
                    {"diff", {1, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p(a[0]->toMatrix());
                        return NumMatrix(p.diff());
                    }}},
                    {"ordadd", {2, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix()), p2(a[1]->toMatrix());
                        return NumMatrix(p1.ordinalAdd(p2));
                    }}},
                    {"ordmul", {2, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix()), p2(a[1]->toMatrix());
                        return NumMatrix(p1.ordinalMul(p2));
                    }}},
                    {"ordsub", {2, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix()), p2(a[1]->toMatrix());
                        return NumMatrix(p1.ordinalSub(p2));
                    }}},
                    {"orddiv", {2, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix()), p2(a[1]->toMatrix());
                        auto res = p1.ordinalDiv(p2);
                        Polynom<Field> div = res.first, mod = res.second;
                        unsigned mx = std::max(div.degree(), mod.degree()) + 1;
                        return NumMatrix(div.toMatrix(mx).joinVertical(mod.toMatrix(mx)));
                    }}},
            };
    string s;
    if(expr.empty())
    {
        s = safeGetline(">> ");
        if(!history_length || s != history_get(history_length)->line)
            add_history(s.c_str());
    }
    else
    {
        s = expr;
    }
    auto v = splitExpression(s);
    static map<char, NumMatrix> mmap;
    vector<pair<token_type, string> > opst;
    vector<NumMatrix> st;
    vector<ssize_t> st_height;
    ssize_t st_size = 0;
    bool dollar = false;
    for(pair<token_type, string> &i : v)
    {
        switch(i.first)
        {
            case TOKEN_POLY:
                ++st_size;
                break;
            case TOKEN_DOLLAR:
                dollar = true;
                break;
            case TOKEN_MATRIX:
                if(dollar)
                {
                    if(i.second[0] == '_')
                        die("Invalid expression");
                    mmap.erase(i.second[0]);
                }
            case TOKEN_NUMBER:
                ++st_size;
                break;
            case TOKEN_OP:
                while(opst.size() && opst.back().first == TOKEN_OP &&
                      priority[int(i.second[0])] + rightassoc[int(i.second[0])] <=
                      priority[int(opst.back().second[0])])
                {
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
                }
                if(st_size < 0)
                    die("Invalid expression");
            case TOKEN_FUNC:
                if(!operations.count(i.second))
                    die("Invalid function: " + i.second);
                if(i.first == TOKEN_FUNC)
                    st_height.push_back(st_size);
            case TOKEN_LEFTPAR:
                opst.push_back(i);
                break;
            case TOKEN_RIGHTPAR:
                while(opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
                }
                if(opst.empty() || st_size <= 0)
                    die("Invalid expression");
                opst.pop_back();
                if(opst.size() && opst.back().first == TOKEN_FUNC)
                {
                    if(operations[opst.back().second].first > 0 &&
                       st_size - operations[opst.back().second].first != st_height.back())
                    {
                        die("Invalid number of arguments: " + opst.back().second);
                    }
                    st_size = st_height.back() + 1;
                    opst.pop_back();
                    st_height.pop_back();
                }
                break;
            case TOKEN_COMMA:
                while(opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
                }
                if(opst.empty() || st_size <= 0)
                    die("Invalid expression");
                break;
        }
        if(i.first != TOKEN_DOLLAR)
            dollar = false;
    }
    while(opst.size())
    {
        if(opst.back().first == TOKEN_LEFTPAR || opst.back().first == TOKEN_RIGHTPAR)
            die("Invalid expression");
        st_size -= operations[opst.back().second].first - 1;
        opst.pop_back();
    }
    if(st_size != 1)
        die("Invalid expression");
    Rational tt;
    expr = "";
    st_height.clear();
    try
    {
        bool matrix_read = false;
        for(pair<token_type, string> &i : v)
        {
            istringstream is, iis;
            switch(i.first)
            {
                case TOKEN_POLY:
                    st.push_back(NumMatrix(Polynom<Field>(i.second)));
                    break;
                case TOKEN_NUMBER:
                    is.str(i.second);
                    is >> tt;
                    if(tt == BigInteger(tt))
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
                    if(!mmap.count(i.second[0]))
                    {
                        if(i.second[0] == '_')
                        {
                            die("_ is not defined");
                        }
                        matrix_read = true;
                        mmap[i.second[0]] = getMatrix<Field>(string("Matrix ") + i.second + ':');
                    }
                    st.push_back(mmap[i.second[0]]);
                    break;
                case TOKEN_OP:
                    while(opst.size() && opst.back().first == TOKEN_OP &&
                          priority[int(i.second[0])] + rightassoc[int(i.second[0])] <=
                          priority[int(opst.back().second[0])])
                    {
                        processOp(opst.back().second, st, 0, operations);
                        opst.pop_back();
                    }
                case TOKEN_FUNC:
                    if(i.first == TOKEN_FUNC)
                        st_height.push_back(st.size());
                case TOKEN_LEFTPAR:
                    opst.push_back(i);
                    break;
                case TOKEN_RIGHTPAR:
                    while(opst.size() && opst.back().first != TOKEN_LEFTPAR)
                    {
                        processOp(opst.back().second, st, 0, operations);
                        opst.pop_back();
                    }
                    opst.pop_back();
                    if(opst.size() && opst.back().first == TOKEN_FUNC)
                    {
                        processOp(opst.back().second, st, st.size() - st_height.back(), operations);
                        st_height.pop_back();
                        opst.pop_back();
                    }
                    break;
                case TOKEN_COMMA:
                    while(opst.size() && opst.back().first != TOKEN_LEFTPAR)
                    {
                        processOp(opst.back().second, st, 0, operations);
                        opst.pop_back();
                    }
                    break;
                case TOKEN_DOLLAR:;
            }
        }
        while(opst.size())
        {
            processOp(opst.back().second, st, 0, operations);
            opst.pop_back();
        }
        mmap['_'] = st[0].toMatrix();
        if((v.size() == 1 || (v.size() == 2 && v[0].first == TOKEN_DOLLAR)) && matrix_read)
            return;
        auto res = st[0].toMatrix();
        cout << res;
        if(res.width() == 1 && res.height() == 1)
            printDecimalResult(res[0][0]);
    }
    catch(matrix_error &e)
    {
        cout << "Error: " << e.what() << endl;
    }
    catch(invalid_number_error &e)
    {
        cout << e.what() << endl;
    }
    catch(end_of_input &e) {}
}


int main(int argc, char **argv)
{
    try
    {
        string arg1 = (argc > 1 ? argv[1] : "");
        string arg2 = (argc > 2 ? argv[2] : "");
        if(isdigit(arg2) && !isdigit(arg1))
            swap(arg1, arg2);
        using_history();
        if(isdigit(arg1))
        {
            _FINITE_ORDER = BigInteger(arg1);
            if(_FINITE_ORDER < 2)
                die("Order must be at least 2");
            /*if(!_FINITE_ORDER.isPrime())
                cout << "WARNING: finite field order is not prime\n";*/
            do
            {
                try
                {
                    f_expr<Finite>(arg2);
                }
                catch(invalid_expression &e) {}
            }
            while(arg2.empty());
        }
        else
        {
            do
            {
                try
                {
                    f_expr<Rational>(arg1);
                }
                catch(invalid_expression &e) {}
            }
            while(arg1.empty());
        }
    }
    catch(matrix_error e)
    {
        cout << "Matrix error: " << e.what() << endl;
        return 1;
    }
    catch(end_of_input &e)
    {
        return 0;
    }
    return 0;
}
