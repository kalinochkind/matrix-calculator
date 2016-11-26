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

template<class NumMatrix>
NumMatrix getMatrix(string prompt)
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
    NumMatrix m(height, width);
    istringstream is;
    is.str(sum);
    is >> m;
    return m;
}

template<class NumMatrix>
NumMatrix f_solve(vector<NumMatrix> a)
{
    unsigned sz = a[0].height();
    if (!sz || a[0].width() != sz + 1)
        die("Invalid use of solve: N*N+1 matrix required");
    NumMatrix sys = a[0].submatrix(0, 0, sz - 1, sz - 1);
    NumMatrix right = a[0].submatrix(0, sz, sz - 1, sz);
    sys.inverseExt(right);
    return right.transposed();

}

template<class NumMatrix>
void processOp(string op, vector<NumMatrix> &st, map<string, pair<int, NumMatrix (*)(vector<NumMatrix>)> > &operations)
{
    if (operations[op].first == 1)
    {
        st.back() = operations[op].second({st.back()});
    }
    else if (operations[op].first == 2)
    {
        NumMatrix a = st.back();
        st.pop_back();
        st.back() = operations[op].second({st.back(), a});
    }
    else
    {
        NumMatrix a = st.back();
        st.pop_back();
        NumMatrix b = st.back();
        st.pop_back();
        st.back() = operations[op].second({st.back(), b, a});
    }
}

template<class Field>
void f_expr()
{

    typedef Matrix<Field> NumMatrix;
    map<string, pair<int, NumMatrix (*)(vector<NumMatrix>)> > operations =
            {
                    {"+",      {2, [](vector<NumMatrix> a) { return a[0] + a[1]; }}},
                    {"^",      {2, [](vector<NumMatrix> a) {
                        if (a[1].width() != 1 || a[1].height() != 1 || a[1][0][0] != int(a[1][0][0]))
                            die("Invalid use of ^: integer required");
                        return a[0].power(int(a[1][0][0]));
                    }}},
                    {"*",      {2, [](vector<NumMatrix> a) {
                        if (a[0].height() == 1 && a[0].width() == 1)
                            return a[1] * a[0][0][0];
                        else if (a[1].height() == 1 && a[1].width() == 1)
                            return a[0] * a[1][0][0];
                        else
                            return a[0] * a[1];
                    }}},
                    {"/",      {2, [](vector<NumMatrix> a) {
                        if (a[0].height() == 1 && a[0].width() == 1)
                            return a[1].inverted() * a[0][0][0];
                        else if (a[1].height() == 1 && a[1].width() == 1)
                            return a[0] * a[1].inverted()[0][0];
                        else
                            return a[0] * a[1].inverted();
                    }}},
                    {"-",      {2, [](vector<NumMatrix> a) { return a[0] - a[1]; }}},
                    {"_",      {1, [](vector<NumMatrix> a) { return -a[0]; }}},
                    {"det",    {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].det()); }}},
                    {"rank",   {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].rank()); }}},
                    {"trace",  {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].trace()); }}},
                    {"t",      {1, [](vector<NumMatrix> a) { return a[0].transposed(); }}},
                    {"inv",    {1, [](vector<NumMatrix> a) { return a[0].inverted(); }}},
                    {"id",     {1, [](vector<NumMatrix> a) {
                        if (a[0].width() != 1 || a[0].height() == 1)
                            die("Invalid use of id");
                        return NumMatrix::identity(abs(int(a[0][0][0])));
                    }}},
                    {"=",      {2, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0] == a[1]); }}},
                    {"width",  {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].width()); }}},
                    {"height", {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].height()); }}},
                    {"solve",  {1, f_solve}},
                    {"at",     {3, [](vector<NumMatrix> a) {
                        if (a[1].width() != 1 || a[1].height() != 1 || a[2].width() != 1 || a[2].height() != 1)
                            die("Invalid use of at");
                        if (int(a[1][0][0]) != a[1][0][0] || int(a[2][0][0]) != a[2][0][0])
                            die("at: indices must be integers");
                        if (int(a[1][0][0]) < 0 || int(a[1][0][0]) >= int(a[0].height()) ||
                            int(a[2][0][0]) < 0 || int(a[2][0][0]) >= int(a[0].width()))
                            die("at: out of range");
                        return NumMatrix::fromNumber(a[0][int(a[1][0][0])][int(a[2][0][0])]);
                    }}}
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
                if (st_size <= 0)
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
        NumMatrix tt(1, 1);
        istringstream is;
        switch (i.first)
        {
            case TOKEN_NUMBER:
                is.str(i.second);
                is >> tt[0][0];
                st.push_back(tt);
                break;
            case TOKEN_MATRIX:
                if (!mmap.count(i.second[0]))
                    mmap[i.second[0]] = getMatrix<NumMatrix>(string("Matrix ") + i.second + ':');
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
    cout << "Result:\n" << st[0];
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
            if(_FINITE_ORDER < 2)
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
