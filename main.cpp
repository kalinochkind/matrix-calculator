#include "matrix.h"
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

pair<unsigned, unsigned> getSize()
{
    string s;
    cout << "Size: ";
    while (s.empty() || s[0] == ' ' || s[0] == '\n')
        getline(cin, s);
    unsigned a, b;
    istringstream iss(s);
    if (!(iss >> a))
        return {0, 0};
    if (!(iss >> b))
        return {a, a};
    return {a, b};
};

NumMatrix getMatrix(unsigned szx = 0, unsigned szy = 0)
{
    pair<unsigned, unsigned> sz;
    if (szx > 0 && szy > 0)
    {
        sz.first = szx;
        sz.second = szy;
    }
    else
        sz = getSize();
    NumMatrix m(sz.first, sz.second);
    cin >> m;
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
    NumMatrix a = getMatrix();
    NumMatrix b = getMatrix();
    cout << "Product:\n" << (a * b);
}

void f_solve()
{
    unsigned sz = getSize().first;
    if (!sz)
        return;
    NumMatrix a = getMatrix(sz, sz + 1);
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

map<string, void (*)()> ops = {{"rank",  f_rank},
                               {"det",   f_det},
                               {"inv",   f_inv},
                               {"mul",   f_mul},
                               {"solve", f_solve}};

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        cout << "Usage: matrix OPERATION\n\n";
        cout << "Operations: ";
        int t = 0;
        for (auto &i : ops)
        {
            cout << (t++ ? ", " : "") << i.first;
        }
        cout << endl << endl;
        return 0;
    }
    if (ops[argv[1]])
    {
        ops[argv[1]]();
    }
    else
    {
        cout << "Invalid operation" << endl;
        return 1;
    }

    return 0;
}