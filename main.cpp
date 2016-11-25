#include "matrix.h"
#include "parser.h"
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

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
        {
            cout << "Incorrect matrix\n";
            exit(1);
        }
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

NumMatrix f_solve(vector<NumMatrix> a)
{
    unsigned sz = a[0].height();
    if (!sz || a[0].width() != sz + 1)
    {
        cout << "Invalid use of solve: N*N+1 matrix required\n";
        exit(1);
    }
    NumMatrix sys = a[0].submatrix(0, 0, sz - 1, sz - 1);
    NumMatrix right = a[0].submatrix(0, sz, sz - 1, sz);
    sys.inverseExt(right);
    return right.transposed();

}

map<string, pair<int, NumMatrix (*)(vector<NumMatrix>)> > operations =
        {
                {"+",      {2, [](vector<NumMatrix> a) { return a[0] + a[1]; }}},
                {"^",      {2, [](vector<NumMatrix> a) {
                    if(a[1].width() != 1 || a[1].height() != 1 || a[1][0][0].denominator() != 1)
                    {
                        cout << "Invalid use of ^: integer required\n";
                        exit(0);
                    }
                    return a[0].power(int(a[1][0][0].numerator()));
                }}},
                {"*",      {2, [](vector<NumMatrix>a) {
                    if (a[0].height() == 1 && a[0].width() == 1)
                        return a[1] * a[0][0][0];
                    else if (a[1].height() == 1 && a[1].width() == 1)
                        return a[0] * a[1][0][0];
                    else
                        return a[0] * a[1];
                }}},
                {"-",      {2, [](vector<NumMatrix> a) { return a[0] - a[1]; }}},
                {"_",      {1, [](vector<NumMatrix> a) { return -a[0]; }}},
                {"det",    {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].det()); }}},
                {"rank",   {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].rank()); }}},
                {"trace",  {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].trace()); }}},
                {"t",      {1, [](vector<NumMatrix> a) { return a[0].transposed(); }}},
                {"inv",    {1, [](vector<NumMatrix> a) { return a[0].inverted(); }}},
                {"id",     {1, [](vector<NumMatrix> a) {
                    assert(a[0].width() == 1 && a[0].height() == 1);
                    return NumMatrix::identity(abs(int(a[0][0][0].numerator())));
                }}},
                {"=",      {2, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0] == a[1]); }}},
                {"width",  {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].width()); }}},
                {"height", {1, [](vector<NumMatrix> a) { return NumMatrix::fromNumber(a[0].height()); }}},
                {"solve", {1, f_solve}},
                {"at", {3, [](vector<NumMatrix> a) {
                        if(a[1].width() != 1 || a[1].height() != 1 || a[2].width() != 1 || a[2].height() != 1)
                        {
                            cout << "Invalid use of at\n";
                            exit(1);
                        }
                        if(a[1][0][0].denominator() != 1 || a[2][0][0].denominator() != 1)
                        {
                            cout << "at: indices must be integers\n";
                            exit(1);
                        }
                        if(a[1][0][0] < 0 || a[1][0][0] >= a[0].height() || a[2][0][0] < 0 || a[2][0][0] >= a[0].width())
                        {
                            cout << "at: out of range\n";
                            exit(1);
                        }
                        return NumMatrix::fromNumber(a[0][int(a[1][0][0].numerator())][int(a[2][0][0].numerator())]);
                    }}}
        };

void processOp(string op, vector<NumMatrix> &st)
{
    if (operations[op].first == 1)
    {
        st.back() = operations[op].second({st.back()});
    }
    else if(operations[op].first == 2)
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

void f_expr()
{
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
                if(st_size <= 0)
                {
                    std::cout << "Invalid expression" << endl;
                    exit(1);
                }
            case TOKEN_FUNC:
                if(!operations.count(i.second))
                {
                    cout << "Invalid function: " << i.second << endl;
                    exit(1);
                }
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
                {
                    cout << "Invalid expression" << endl;
                    exit(1);
                }
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
                    processOp(opst.back().second, st);
                    opst.pop_back();
                }
                if (opst.empty() || st_size <= 0)
                {
                    cout << "Invalid expression" << endl;
                    exit(1);
                }
                break;
        }
    }
    while (opst.size())
    {
        if (opst.back().first == TOKEN_LEFTPAR || opst.back().first == TOKEN_RIGHTPAR)
        {
            cout << "Invalid expression" << endl;
            exit(1);
        }
        st_size -= operations[opst.back().second].first - 1;
        opst.pop_back();
    }
    if(st_size != 1)
    {
        cout << "Invalid expression" << endl;
        exit(1);
    }
    for (pair<token_type, string> &i : v)
    {
        NumMatrix tt(1, 1);
        switch (i.first)
        {
            case TOKEN_NUMBER:
                tt[0][0] = Rational(i.second);
                st.push_back(tt);
                break;
            case TOKEN_MATRIX:
                if (!mmap.count(i.second[0]))
                    mmap[i.second[0]] = getMatrix(string("Matrix ") + i.second + ':');
                st.push_back(mmap[i.second[0]]);
                break;
            case TOKEN_OP:
                while (opst.size() && opst.back().first == TOKEN_OP &&
                       priority[int(i.second[0])] + rightassoc[int(i.second[0])] <=
                       priority[int(opst.back().second[0])])
                {
                    processOp(opst.back().second, st);
                    opst.pop_back();
                }
            case TOKEN_FUNC:
            case TOKEN_LEFTPAR:
                opst.push_back(i);
                break;
            case TOKEN_RIGHTPAR:
                while (opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    processOp(opst.back().second, st);
                    opst.pop_back();
                }
                opst.pop_back();
                if (opst.size() && opst.back().first == TOKEN_FUNC)
                {
                    processOp(opst.back().second, st);
                    opst.pop_back();
                }
                break;
            case TOKEN_COMMA:
                while (opst.size() && opst.back().first != TOKEN_LEFTPAR)
                {
                    processOp(opst.back().second, st);
                    opst.pop_back();
                }
                break;
        }
    }
    while (opst.size())
    {
        processOp(opst.back().second, st);
        opst.pop_back();
    }
    cout << "Result:\n" << st[0];
}


int main()
{

    f_expr();
    return 0;
}
