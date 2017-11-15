#include "matrix.h"
#include "finite.h"
#include "complex.h"
#include <sstream>

const unsigned strassen_threshold = 32;

template<class Field>
void Matrix<Field>::_strassen(unsigned N, Field *a, Field *b, Field *res)
{
    if(N <= strassen_threshold)
    {
        for(unsigned i = 0; i < N; ++i)
        {
            for(unsigned j = 0; j < N; ++j)
            {
                for(unsigned q = 0; q < N; ++q)
                {
                    res[i * N + j] += a[i * N + q] * b[q * N + j];
                }
            }
        }
        return;
    }
    Field *buf = new Field[N * N / 2];
    Field *p = new Field[N * N / 4 * 7];
#define AQUAD(x, y) a[(x?N*N/2:0)+(y?N/2:0)+i*N+j]
#define BQUAD(x, y) b[(x?N*N/2:0)+(y?N/2:0)+i*N+j]
#define P(n) p[N*N/4*n+i*N/2+j]
#define FORFOR for (unsigned i = 0; i < N / 2; ++i) for (unsigned j = 0; j < N / 2; ++j)
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 0) + AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 0) + BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(1, 0) + AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 0);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 0);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 1) - BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 2 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(1, 0) - BQUAD(0, 0);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 3 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 0) + AQUAD(0, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 4 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(1, 0) - AQUAD(0, 0);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(0, 0) + BQUAD(0, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 5 * N * N / 4);
    FORFOR
        {
            buf[i * N / 2 + j] = AQUAD(0, 1) - AQUAD(1, 1);
            buf[N * N / 4 + i * N / 2 + j] = BQUAD(1, 0) + BQUAD(1, 1);
        }
    _strassen(N / 2, buf, buf + N * N / 4, p + 6 * N * N / 4);
    FORFOR
        {
            res[N * i + j] = P(0) + P(3) - P(4) + P(6);
            res[N * i + N / 2 + j] = P(2) + P(4);
            res[N * N / 2 + N * i + j] = P(1) + P(3);
            res[N * N / 2 + N * i + N / 2 + j] = P(0) - P(1) + P(2) + P(5);
        }
    delete[] buf;
    delete[] p;
#undef AQUAD
#undef BQUAD
#undef P
#undef FORFOR
}


template<class Field>
const Matrix<Field> Matrix<Field>::operator*(const Matrix &a) const
{
    if(N != a.M)
        throw matrix_error("Trying to multiply matrices of different size");
    if(std::max(M, a.N) < strassen_threshold)
    {
        return _multiplyCube(a);
    }
    else
    {
        return _multiplyStrassen(a);
    }
}

template <class Field>
Matrix<Field> Matrix<Field>::_multiplyCube(const Matrix &a) const
{
    assert(N == a.M);
    Matrix res(M, a.N);
    for(unsigned i = 0; i < M; ++i)
    {
        for(unsigned j = 0; j < a.N; ++j)
        {
            for(unsigned q = 0; q < N; ++q)
            {
                res[i][j] += (*this)[i][q] * a[q][j];
            }
        }
    }
    return res;
}

template<class Field>
const Matrix<Field> Matrix<Field>::_multiplyStrassen(const Matrix &m) const
{
    assert(N == m.M);
    Matrix res(M, m.N);
    unsigned mx = std::max(M, std::max(N, m.N));
    while(mx & (mx - 1))
        ++mx;
    Field *a = new Field[mx * mx];
    Field *b = new Field[mx * mx];
    Field *r = new Field[mx * mx];
    for(unsigned i = 0; i < M; ++i)
        for(unsigned j = 0; j < N; ++j)
            a[i * mx + j] = (*this)[i][j];
    for(unsigned i = 0; i < N; ++i)
        for(unsigned j = 0; j < m.N; ++j)
            b[i * mx + j] = m[i][j];
    _strassen(mx, a, b, r);
    for(unsigned i = 0; i < M; ++i)
        for(unsigned j = 0; j < m.N; ++j)
            res[i][j] = r[i * mx + j];
    delete[] a;
    delete[] b;
    delete[] r;
    return res;
}

template<class Field>
Matrix<Field>::Matrix(const Matrix &a): Matrix(a.M, a.N)
{
    for(unsigned i = 0; i < M * N; ++i)
    {
        arr[i] = a.arr[i];
    }
}

template<class Field>
const Matrix<Field> Matrix<Field>::fromRow(const std::string &s)
{
    std::istringstream iss;
    iss.str(s);
    Field f;
    std::vector<Field> v;
    while(iss >> f)
        v.push_back(f);
    if(v.empty())
        return Matrix(1, 1);
    Matrix m(1, v.size());
    for(unsigned i = 0; i < v.size(); ++i)
        m[0][i] = v[i];
    return m;
}

template<class Field>
Matrix<Field> &Matrix<Field>::operator=(const Matrix &a)
{
    if(arr == a.arr)
        return *this;
    M = a.M;
    N = a.N;
    delete[] arr;
    arr = new Field[M * N];
    for(unsigned i = 0; i < M * N; ++i)
    {
        arr[i] = a.arr[i];
    }
    return *this;
}

template<class Field>
unsigned Matrix<Field>::gauss(Matrix *ext)
{
    if(ext)
    {
        assert(M == ext->M);
    }
    unsigned row = 0;
    for(unsigned col = 0; col < N; ++col)
    {
        for(unsigned i = row; i < M; ++i)
        {
            if((*this)[i][col])
            {
                if(i != row)
                {
                    for(unsigned j = 0; j < N; ++j)
                    {
                        std::swap((*this)[i][j], (*this)[row][j]);
                        (*this)[i][j] = -(*this)[i][j];
                    }
                    if(ext)
                    {
                        for(unsigned j = 0; j < ext->N; ++j)
                        {
                            std::swap((*ext)[i][j], (*ext)[row][j]);
                            (*ext)[i][j] = -(*ext)[i][j];
                        }
                    }
                }
                break;
            }
        }
        if(!(*this)[row][col])
        {
            continue;
        }
        for(unsigned i = row + 1; i < M; ++i)
        {
            Field coeff = (*this)[i][col] / (*this)[row][col];
            for(unsigned j = 0; j < N; ++j)
            {
                (*this)[i][j] -= (*this)[row][j] * coeff;
            }
            if(ext)
            {
                for(unsigned j = 0; j < ext->N; ++j)
                    (*ext)[i][j] -= (*ext)[row][j] * coeff;
            }
        }
        if(++row == M)
            return row;
    }
    return row;
}

template<class Field>
const Matrix<Field> Matrix<Field>::identity(unsigned n)
{
    Matrix res(n);
    for(unsigned i = 0; i < n; ++i)
    {
        res[i][i] = 1;
    }
    return res;
}

template<class Field>
Matrix<Field> &Matrix<Field>::operator+=(const Matrix &a)
{
    if(M != a.M || N != a.N)
        throw matrix_error("Trying to add matrices of different size");
    for(unsigned i = 0; i < N * M; ++i)
    {
        arr[i] += a.arr[i];
    }
    return *this;
}

template<class Field>
const Matrix<Field> Matrix<Field>::operator-() const
{
    Matrix temp(*this);
    for(unsigned i = 0; i < N * M; ++i)
    {
        temp.arr[i] = -temp.arr[i];
    }
    return temp;
}

template <class Field>
Matrix<Field> &Matrix<Field>::operator-=(const Matrix &a)
{
    if(M != a.M || N != a.N)
        throw matrix_error("Trying to subtract matrices of different size");
    for(unsigned i = 0; i < N * M; ++i)
    {
        arr[i] -= a.arr[i];
    }
    return *this;
}

template <class Field>
Matrix<Field> &Matrix<Field>::operator*=(Field a)
{
    for(unsigned i = 0; i < N * M; ++i)
    {
        arr[i] *= a;
    }
    return *this;
}

template <class Field>
const Field Matrix<Field>::det() const
{
    if(N != M)
        throw matrix_error("Trying to calculate determinant of a non-square matrix");
    Matrix tmp(*this);
    if(tmp.gauss() != N)
        return Field();
    Field ans = 1;
    for(unsigned i = 0; i < N * N; i += N + 1)
    {
        ans *= tmp.arr[i];
    }
    return ans;
}

template <class Field>
const Matrix<Field> Matrix<Field>::transposed() const
{
    Matrix res(N, M);
    for(unsigned i = 0; i < M; ++i)
        for(unsigned j = 0; j < N; ++j)
        {
            res[j][i] = (*this)[i][j];
        }
    return res;
}

template <class Field>
const Field Matrix<Field>::trace() const
{
    if(N != M)
        throw matrix_error("Trying to calculate trace of a non-square matrix");
    Field ans;
    for(unsigned i = 0; i < N * N; i += N + 1)
    {
        ans += arr[i];
    }
    return ans;
}

template <class Field>
const Matrix<Field> Matrix<Field>::inverted() const
{
    Matrix ext = Matrix::identity(N);
    Matrix(*this).inverseExt(ext);
    return ext;
}

template <class Field>
bool Matrix<Field>::operator==(const Matrix &a) const
{
    if(M != a.M || N != a.N)
        return false;
    for(unsigned i = 0; i < N * M; ++i)
    {
        if(arr[i] != a.arr[i])
            return false;
    }
    return true;
}

template <class Field>
const Matrix<Field> Matrix<Field>::submatrix(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const
{
    if(x1 > x2)
        std::swap(x1, x2);
    if(y1 > y2)
        std::swap(y1, y2);
    Matrix res(x2 - x1 + 1, y2 - y1 + 1);
    for(unsigned i = x1; i <= x2; ++i)
    {
        for(unsigned j = y1; j <= y2; ++j)
        {
            res[i - x1][j - y1] = (*this)[i][j];
        }
    }
    return res;
}

template <class Field>
void Matrix<Field>::inverseExt(Matrix &ext)
{
    if(N != M)
        throw matrix_error("Trying to inverse a non-square matrix");
    if(M != ext.M)
    {
        throw matrix_error("Invalid use of inverseExt");
    }
    if(gauss(&ext) != N)
    {
        throw matrix_error("Cannot inverse a singular matrix");
    }
    for(unsigned i = N - 1; i != unsigned(-1); --i)
    {
        Field coeff = Field(1) / (*this)[i][i];
        for(unsigned j = 0; j < N; ++j)
        {
            (*this)[i][j] *= coeff;
        }
        for(unsigned j = 0; j < ext.N; ++j)
        {
            ext[i][j] *= coeff;
        }
        for(unsigned j = i + 1; j < N; j++)
        {
            if(!(*this)[i][j])
                continue;
            for(unsigned q = 0; q < ext.N; q++)
            {
                ext[i][q] -= ext[j][q] * (*this)[i][j];
            }
            (*this)[i][j] = Field();
        }
    }
}

template <class Field>
const Matrix<Field> Matrix<Field>::fundamental() const
{
    Matrix m(*this);
    unsigned rk = m.gauss();
    std::vector<char> dependent(N);
    std::vector<int> dep_height(N);
    int last_nonzero = -1;
    for(unsigned i=0;i<N;++i)
    {
        for(int j=M-1;j>last_nonzero;--j)
        {
            if(m[j][i])
            {
                last_nonzero = j;
                dependent[i] = 1;
                dep_height[i] = j;
                break;
            }
        }
    }
    for(unsigned i=0;i<N;++i)
    {
        if(!dependent[i])
            continue;
        Field inv = Field(1) / m[dep_height[i]][i];
        for(unsigned j=i;j<N;++j)
        {
            m[dep_height[i]][j] *= inv;
        }
        for(int j=0;j<dep_height[i];++j)
        {
            if(!m[j][i])
                continue;
            for(unsigned q=i+1;q<N;++q)
            {
                m[j][q] -= m[j][i] * m[dep_height[i]][q];
            }
            m[j][i] = 0;
        }
    }
    Matrix ans(N, N - rk);
    int cur = 0;
    unsigned col = -1;
    for(unsigned i=0;i<N;++i)
    {
        if(dependent[i])
            continue;
        do
        {
            ++col;
        }
        while(dependent[col]);
        for(unsigned j=0;j<N;++j)
        {
            if(dependent[j])
            {
                ans[j][cur] = -m[dep_height[j]][i];
            }
            else if(j == col)
            {
                ans[j][cur] = 1;
            }
        }
        ++cur;
    }
    makeIntColumns(ans);
    return ans;
}

template <class Field>
const Matrix<Field> Matrix<Field>::partial() const
{
    if(N < 2)
        throw matrix_error("Invalid matrix");
    Matrix sys = submatrix(0, 0, M - 1, N - 2);
    Matrix right = submatrix(0, N - 1, M - 1, N - 1);
    Matrix ans(1, N - 1);
    if(sys.gauss(&right) != rank())
    {
        throw matrix_error("No solutions");
    }
    int row = M - 1;
    for(int col=N-2;col>=0;--col)
    {
        while(row >= 0 && !sys[row][col])
            --row;
        if(row < 0)
            return ans;
        Field cval = right[row][0] / sys[row][col];
        if(!cval)
            continue;
        ans[0][col] = cval;
        for(int i=0;i<=row;++i)
        {
            right[i][0] -= cval * sys[i][col];
        }
    }
    return ans;

}

template <class Field>
const Matrix<Field> Matrix<Field>::power(const BigInteger &pow) const
{
    if(M != N)
        throw matrix_error("Power is only defined for square matrices");
    Matrix t = identity(width());
    bool t_set = false;
    Matrix a(*this);
    BigInteger p = abs(pow);
    while(p)
    {
        if(p.odd())
        {
            if (t_set)
                t *= a;
            else
                t = a;
            t_set = true;
        }
        p._divide_by_2();
        if(p)
            a *= a;
    }
    return pow < 0 ? t.inverted() : t;
}

template <class Field>
const Matrix<Field> Matrix<Field>::joinHorizontal(const Matrix &a) const
{
    if(M != a.M)
    {
        throw matrix_error("Trying to join matrices of different height");
    }
    Matrix res(M, N + a.N);
    for(unsigned i = 0; i < M; i++)
    {
        for(unsigned j = 0; j < N; j++)
        {
            res[i][j] = (*this)[i][j];
        }
        for(unsigned j = 0; j < a.N; j++)
        {
            res[i][N + j] = a[i][j];
        }
    }
    return res;
}

template <class Field>
const Matrix<Field> Matrix<Field>::joinVertical(const Matrix &a) const
{
    if(N != a.N)
    {
        throw matrix_error("Trying to join matrices of different width");
    }
    Matrix res(M + a.M, N);
    for(unsigned i = 0; i < M; i++)
    {
        for(unsigned j = 0; j < N; j++)
        {
            res[i][j] = (*this)[i][j];
        }
    }
    for(unsigned i = 0; i < a.M; i++)
    {
        for(unsigned j = 0; j < N; j++)
        {
            res[M + i][j] = a[i][j];
        }
    }
    return res;
}

template<class Field>
const Matrix<Field> Matrix<Field>::charPolynom() const
{
    if(N != M)
    {
        throw matrix_error("Characteristic polynom is only defined for square matrices");
    }
    Matrix A(*this), B(identity(N));
    Matrix ans(1, N + 1);
    ans[0][0] = 1;
    for(unsigned i=1;i<=N;++i)
    {
        A = *this * B;
        Field p = A.trace() / i;
        B = A - identity(N) * p;
        ans[0][i] = -p;
    }
    return ans;
}

template<class Field>
void Matrix<Field>::_write_to_ostream(std::ostream &out) const
{
    if(!N || !M)
        return;
    std::vector<std::vector<std::string>> strs(M);
    std::vector<size_t> max_width(N);
    for(unsigned i = 0; i < M; ++i)
    {
        for(unsigned j = 0; j < N; ++j)
        {
            std::string s = toString((*this)[i][j]);
            strs[i].push_back(s);
            max_width[j] = std::max(max_width[j], s.length());
        }
    }
    for(unsigned i = 0; i < M; ++i)
    {
        for(unsigned j = 0; j < N; ++j)
        {
            if(j)
                out << ' ';
            for(unsigned t=0;t+strs[i][j].length()<max_width[j];++t)
            {
                out << ' ';
            }
            out << strs[i][j];
        }
        out << '\n';
    }
}

template<class Field>
void Matrix<Field>::_read_from_istream(std::istream &in)
{
    for(unsigned i = 0; i < M; ++i)
    {
        for(unsigned j = 0; j < N; ++j)
        {
            in >> (*this)[i][j];
        }
    }
}

void makeIntColumns(Matrix<Finite> &) {}

void makeIntColumns(Matrix<Rational> &m)
{
    for(unsigned col=0;col<m.width();++col)
    {
        BigInteger g = 1;
        unsigned pos = 0, neg = 0;
        for(unsigned i=0;i<m.height();++i)
        {
            if(m[i][col].denominator() != 1)
                g *= m[i][col].denominator() / gcd(g, m[i][col].denominator());
            pos += m[i][col] > 0;
            neg += m[i][col] < 0;
        }
        if(neg > pos)
            g = -g;
        for(unsigned i=0;i<m.height();++i)
        {
            m[i][col] *= g;
        }
    }
}

void makeIntColumns(Matrix<Complex> &m)
{
    for(unsigned col=0;col<m.width();++col)
    {
        BigInteger g = 1;
        unsigned pos = 0, neg = 0;
        for(unsigned i=0;i<m.height();++i)
        {
            if(m[i][col].re().denominator() != 1 || m[i][col].im().denominator() != 1 )
            {
                g *= m[i][col].denominator() / gcd(g, m[i][col].denominator());
            }
            pos += m[i][col] > 0;
            neg += m[i][col] < 0;
        }
        if(neg > pos)
            g = -g;
        for(unsigned i=0;i<m.height();++i)
        {
            m[i][col] *= g;
        }
    }
}


template class Matrix<Rational>;
template class Matrix<Finite>;
template class Matrix<Complex>;
