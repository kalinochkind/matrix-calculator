#include "matrix.h"
#include "parser.h"
#include "finite.h"
#include "polynom.h"
#include <map>
#include <set>
#include <sstream>
#include <cstdlib>

using namespace std;

void die(const string &s)
{
    cout << s << endl;
    exit(1);
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

string safeGetline()
{
    string s;
    while(s.empty())
    {
        getline(cin, s);
        if(cin.eof())
        {
            cout << "\n";
            exit(0);
        }
        while(s.size() && isspace(s.back()))
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
    while(s.length())
    {
        istringstream is;
        is.str(s);
        unsigned cwidth = 0;
        Rational dummy;
        while(is >> dummy)
            ++cwidth;
        if(!cwidth)
            break;
        if(width && width != cwidth)
            throw matrix_error("Incorrect matrix");
        width = cwidth;
        ++height;
        sum += s + ' ';
        getline(cin, s);
        if(cin.eof())
            exit(0);
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
    if(operations[op].first == 1)
    {
        st.back() = operations[op].second({&st.back()});
    }
    else if(operations[op].first == 2)
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
        return toMatrix() + m.toMatrix();
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

Matrix<Finite> f_cfrac(Finite a)
{
    return Matrix<Finite>::fromNumber(a);
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


template<class T>
T f_revcfrac(const Matrix<T> &a)
{
    T res = a[0][a.width() - 1];
    for(int i = int(a.width()) - 2; i >= 0; --i)
    {
        res = T(1) / res + a[0][i];
    }
    return res;
}

void printDecimalResult(const Finite &) {}

void printDecimalResult(const Rational &a)
{
    if(a.denominator() == 1)
        return;
    cout << "\nDecimal:\n" << a.asDecimal(30) << endl;
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
                        Matrix<Field> m = a[0]->toMatrix();
                        unsigned sz = m.width() - 1;
                        //if (!sz || a[0]->toMatrix().width() != sz + 1)
                        //    die("Invalid use of solve: N*N+1 matrix required");
                        Matrix<Field> sys = m.submatrix(0, 0, m.height() - 1, sz - 1);
                        Matrix<Field> right = m.submatrix(0, sz, m.height() - 1, sz);
                        if(m.rank() != sys.rank())
                        {
                            throw matrix_error("No solutions");
                        }
                        sys.inverseExt(right);
                        return NumMatrix(right.transposed());
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
                            die("revcfrac: matrix 1*n required");
                        return NumMatrix(f_revcfrac(m));
                    }}},
                    {"joinh",     {2, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().joinHorizontal(a[1]->toMatrix()));
                    }}},
                    {"joinv",     {2, [](const vector<NumMatrix *> &a) {
                        return NumMatrix(a[0]->toMatrix().joinVertical(a[1]->toMatrix()));
                    }}},
                    {"gauss",     {1, [](const vector<NumMatrix *> &a) {
                        Matrix<Field> m = a[0]->toMatrix();
                        m.gauss();
                        return NumMatrix(m);
                    }}},
                    {"gcd",       {2, [](const vector<NumMatrix *> &a) {
                        Polynom<Field> p1(a[0]->toMatrix()), p2(a[1]->toMatrix());
                        if(p1.degree() < 1 && p2.degree() < 1)
                        {
                            BigInteger i1(p1.toMatrix()[0][0]);
                            BigInteger i2(p2.toMatrix()[0][0]);
                            return NumMatrix(gcd(i1, i2));
                        }
                        return NumMatrix(p1.gcd(p2));
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
            };
    string s;
    if(expr.empty())
    {
        cout << "Expression: ";
        s = safeGetline();
    }
    else
    {
        s = expr;
    }
    auto v = splitExpression(s);
    map<char, NumMatrix> mmap;
    set<char> repeated;
    vector<pair<token_type, string> > opst;
    vector<NumMatrix> st;
    int st_size = 0;
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
                    repeated.insert(i.second[0]);
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
                    st_size -= operations[opst.back().second].first - 1;
                    opst.pop_back();
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
    do
    {
        if(expr.empty())
            cout << endl;
        expr = "";
        try
        {
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
                            mmap[i.second[0]] = getMatrix<Matrix<Field>>(string("Matrix ") + i.second + ':');
                        st.push_back(mmap[i.second[0]]);
                        break;
                    case TOKEN_OP:
                        while(opst.size() && opst.back().first == TOKEN_OP &&
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
                        while(opst.size() && opst.back().first != TOKEN_LEFTPAR)
                        {
                            processOp(opst.back().second, st, operations);
                            opst.pop_back();
                        }
                        opst.pop_back();
                        if(opst.size() && opst.back().first == TOKEN_FUNC)
                        {
                            processOp(opst.back().second, st, operations);
                            opst.pop_back();
                        }
                        break;
                    case TOKEN_COMMA:
                        while(opst.size() && opst.back().first != TOKEN_LEFTPAR)
                        {
                            processOp(opst.back().second, st, operations);
                            opst.pop_back();
                        }
                        break;
                    case TOKEN_DOLLAR:;
                }
            }
            while(opst.size())
            {
                processOp(opst.back().second, st, operations);
                opst.pop_back();
            }
            auto res = st[0].toMatrix();
            cout << "Result:\n" << res;

            if(res.width() == 1 && res.height() == 1)
                printDecimalResult(res[0][0]);
        }
        catch(matrix_error e)
        {
            cout << "Error: " << e.what() << endl;
        }
        st.clear();
        opst.clear();
        for(char i : repeated)
        {
            mmap.erase(i);
        }
    } while(!repeated.empty());

}


int main(int argc, char **argv)
{
    try
    {
        string arg1 = (argc > 1 ? argv[1] : "");
        string arg2 = (argc > 2 ? argv[2] : "");
        if(isdigit(arg2) && !isdigit(arg1))
            swap(arg1, arg2);
        if(isdigit(arg1))
        {
            _FINITE_ORDER = BigInteger(arg1);
            if(_FINITE_ORDER < 2)
                die("Order must be at least 2");
            /*if(!_FINITE_ORDER.isPrime())
                cout << "WARNING: finite field order is not prime\n";*/
            f_expr<Finite>(arg2);
        }
        else
        {
            f_expr<Rational>(arg1);
        }
    }
    catch(matrix_error e)
    {
        cout << "Matrix error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
