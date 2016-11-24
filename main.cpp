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

NumMatrix getMatrix(string prompt = "Enter the matrix:")
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

void f_rank()
{
    NumMatrix a = getMatrix();
    cout << "Rank: " << a.rank() << endl;
}

void f_det()
{
    NumMatrix a = getMatrix();
    cout << "Determinant: " << a.det() << endl;
}

void f_inv()
{
    NumMatrix a = getMatrix();
    cout << "Inverted matrix:\n" << a.inverted();
}

void f_mul()
{
    NumMatrix a = getMatrix("First matrix:");
    NumMatrix b = getMatrix("Second matrix:");
    cout << "Product:\n" << (a * b);
}

void f_solve()
{
    NumMatrix a = getMatrix();
    unsigned sz = a.height();
    if (!sz || a.width() != sz + 1)
    {
        cout << "Invalid matrix\n";
        return;
    }
    NumMatrix sys = a.submatrix(0, 0, sz - 1, sz - 1);
    NumMatrix right = a.submatrix(0, sz, sz - 1, sz);
    sys.inverseExt(right);
    cout << "Solution:";
    for (unsigned i = 0; i < sz; ++i)
    {
        cout << ' ' << right[i][0];
    }
    cout << endl;

}

void f_pow()
{
    NumMatrix a = getMatrix();
    if (a.width() != a.height())
    {
        cout << "This can be done only for square matrices" << endl;
        return;
    }
    cout << "Power: ";
    unsigned p;
    cin >> p;
    NumMatrix t = a.power(p);
    cout << "Result:\n" << t;
}

void processOp(string op, vector<NumMatrix> &st)
{
    if (op == "+")
    {
        assert(st.size() >= 2);
        NumMatrix a = st.back();
        st.pop_back();
        st.back() += a;
    }
    if (op == "*")
    {
        assert(st.size() >= 2);
        NumMatrix a = st.back();
        st.pop_back();
        if (a.width() == 1 && a.height() == 1)
        {
            st.back() *= a[0][0];
        }
        else if (st.back().width() == 1 && st.back().height() == 1)
        {
            st.back() = a * st.back()[0][0];
        }
        else
        {
            st.back() *= a;
        }
    }
    if (op == "^")
    {
        assert(st.size() >= 2);
        NumMatrix a = st.back();
        assert(a.height() == 1);
        assert(a.width() == 1);
        st.pop_back();
        st.back() = st.back().power(int(a[0][0].numerator()));
    }
    if (op == "-")
    {
        assert(st.size() >= 2);
        NumMatrix a = st.back();
        st.pop_back();
        st.back() -= a;
    }
    if (op == "_")
    {
        assert(st.size() >= 1);
        st.back() = -st.back();
    }
    if (op == "det")
    {
        assert(st.size() >= 1);
        NumMatrix a(1);
        a[0][0] = st.back().det();
        st.back() = a;
    }
    if (op == "rank")
    {
        assert(st.size() >= 1);
        NumMatrix a(1);
        a[0][0] = st.back().rank();
        st.back() = a;
    }
    if (op == "trace")
    {
        assert(st.size() >= 1);
        NumMatrix a(1);
        a[0][0] = st.back().trace();
        st.back() = a;
    }
    if (op == "t")
    {
        assert(st.size() >= 1);
        st.back() = st.back().transposed();
    }
    if (op == "inv")
    {
        assert(st.size() >= 1);
        st.back().inverse();
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
    for (auto i : v)
    {
        NumMatrix tt(1, 1);
        switch (i.first)
        {
            case TOKEN_NUMBER:
                tt[0][0] = BigInteger(i.second);
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
                if (opst.empty())
                {
                    cout << "Invalid expression" << endl;
                    exit(1);
                }
                opst.pop_back();
                if (opst.size() && opst.back().first == TOKEN_FUNC)
                {
                    processOp(opst.back().second, st);
                    opst.pop_back();
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
        processOp(opst.back().second, st);
        opst.pop_back();
    }
    if (st.size() == 1)
        cout << "Result:\n" << st[0];
    else
        cout << "Invalid expression\n";
}


map<string, void (*)()> ops = {{"rank",  f_rank},
                               {"det",   f_det},
                               {"inv",   f_inv},
                               {"mul",   f_mul},
                               {"solve", f_solve},
                               {"pow",   f_pow}};

void print_help()
{
    cout << "Usage: matrix [OPERATION]\n\n";
    cout << "Operations: ";
    int t = 0;
    for (auto &i : ops)
    {
        cout << (t++ ? ", " : "") << i.first;
    }
    cout << endl << endl;
}

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        f_expr();
        return 0;
    }
    else if (string(argv[1]) == "help")
    {
        print_help();
        return 0;
    }
    if (ops[argv[1]])
    {
        ops[argv[1]]();
    }
    else
    {
        cout << "Invalid operation" << endl;
        print_help();
        return 1;
    }
    return 0;
}