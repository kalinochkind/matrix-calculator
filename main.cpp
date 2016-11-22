#include "matrix.h"
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

NumMatrix getMatrix(string prompt = "Enter the matrix:")
{
    string s, sum;
    cout << prompt << endl;
    while (s.empty())
        getline(cin, s);
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


void f_expr()
{
    cout << "Expression: ";
    string s;
    getline(cin, s);
    s += ' ';
    map<char, NumMatrix> mmap;
    NumMatrix m, cm;
    string power;
    for (char i : s)
    {
        if (i == '^')
        {
            continue;
        }
        else if (('0' <= i && i <= '9') || i == '-')
        {
            power.push_back(i);
            continue;
        }
        else if(power.length())
        {
            cm = cm.power(atoi(power.c_str()));
            power = "";
        }
        if ('A' <= i && i <= 'Z')  // matrix name
        {
            if(!m.width())
            {
                m = cm;
            }
            else
            {
                m *= cm;
            }
            if (!mmap.count(i))
            {
                mmap[i] = getMatrix(string("Matrix ") + i + ':');
            }
            cm = mmap[i];
        }
    }
    if(!m.width())
    {
        m = cm;
    }
    else
    {
        m *= cm;
    }
    cout << "Result:\n" <<  m;
}

map<string, void (*)()> ops = {{"rank",  f_rank},
                               {"det",   f_det},
                               {"inv",   f_inv},
                               {"mul",   f_mul},
                               {"solve", f_solve},
                               {"pow",   f_pow},
                               {"expr",  f_expr}};

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