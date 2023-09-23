#include "block.H"
#include "blockM.H"

#include "pTransform.H"

namespace Foam
{

namespace briscola
{

// Constructors

template<class Type>
block<Type>::block(const int l, const int m, const int n)
:
    refCount(),
    l_(l),
    m_(m),
    n_(n),
    v_(nullptr),
    T_(eye)
{
    allocate();
}

template<class Type>
block<Type>::block(const int l, const int m, const int n, const zero)
:
    refCount(),
    l_(l),
    m_(m),
    n_(n),
    v_(nullptr),
    T_(eye)
{
    allocate();

    *this = Zero;
}

template<class Type>
block<Type>::block(const int l, const int m, const int n, const Type& s)
:
    refCount(),
    l_(l),
    m_(m),
    n_(n),
    v_(nullptr),
    T_(eye)
{
    allocate();

    *this = s;
}

template<class Type>
block<Type>::block(const labelVector& d)
:
    refCount(),
    l_(d[0]),
    m_(d[1]),
    n_(d[2]),
    v_(nullptr),
    T_(eye)
{
    allocate();
}

template<class Type>
block<Type>::block(const labelVector& d, const zero)
:
    refCount(),
    l_(d[0]),
    m_(d[1]),
    n_(d[2]),
    v_(nullptr),
    T_(eye)
{
    allocate();

    *this = Zero;
}

template<class Type>
block<Type>::block(const labelVector& d, const Type& s)
:
    refCount(),
    l_(d[0]),
    m_(d[1]),
    n_(d[2]),
    v_(nullptr),
    T_(eye)
{
    allocate();

    *this = s;
}

template<class Type>
block<Type>::block
(
    const label l,
    const label m,
    const label n,
    const List<Type>& v
)
:
    refCount(),
    l_(l),
    m_(m),
    n_(n),
    v_(nullptr),
    T_(eye)
{
    if (l_*m_*n_ != v.size())
    {
        FatalErrorInFunction
            << "Cannot create block with size " << l_*m_*n_
            << " from list with size " << v.size() << endl
            << abort(FatalError);
    }
    else
    {
        allocate();

        if (v_)
        {
            if (contiguous<Type>())
            {
                memcpy(v_, v.cdata(), this->byteSize());
            }
            else
            {
                forAllBlockLinear(*this, i)
                {
                    v_[i] = v[i];
                }
            }
        }
    }
}

template<class Type>
block<Type>::block
(
    const label l,
    const label m,
    const label n,
    const Type* v
)
:
    refCount(),
    l_(l),
    m_(m),
    n_(n),
    v_(nullptr),
    T_(eye)
{
    allocate();

    if (v_)
    {
        if (contiguous<Type>())
        {
            memcpy(v_, v, this->byteSize());
        }
        else
        {
            forAllBlockLinear(*this, i)
            {
                v_[i] = v[i];
            }
        }
    }
}

template<class Type>
block<Type>::block(const labelVector& d, const List<Type>& v)
:
    refCount(),
    l_(d[0]),
    m_(d[1]),
    n_(d[2]),
    v_(nullptr),
    T_(eye)
{
    if (l_*m_*n_ != v.size())
    {
        FatalErrorInFunction
            << "Cannot create block with size " << l_*m_*n_
            << " from list with size " << v.size() << endl
            << abort(FatalError);
    }
    else
    {
        allocate();

        if (v_)
        {
            if (contiguous<Type>())
            {
                memcpy(v_, v.cdata(), this->byteSize());
            }
            else
            {
                forAllBlockLinear(*this, i)
                {
                    v_[i] = v[i];
                }
            }
        }
    }
}

template<class Type>
block<Type>::block(const labelVector& d, const Type* v)
:
    refCount(),
    l_(d[0]),
    m_(d[1]),
    n_(d[2]),
    v_(nullptr),
    T_(eye)
{
    allocate();

    if (v_)
    {
        if (contiguous<Type>())
        {
            memcpy(v_, v, this->byteSize());
        }
        else
        {
            forAllBlockLinear(*this, i)
            {
                v_[i] = v[i];
            }
        }
    }
}

template<class Type>
block<Type>::block(const block<Type>& M)
:
    refCount(),
    l_(M.l()),
    m_(M.m()),
    n_(M.n()),
    v_(nullptr),
    T_(M.T())
{
    if (M.v())
    {
        allocate();

        forAllBlockLinear(*this, i)
        {
            v_[i] = M.v()[i];
        }
    }
}

template<class Type>
block<Type>::block(const block<Type>& M, const zero&)
:
    refCount(),
    l_(M.l()),
    m_(M.m()),
    n_(M.n()),
    v_(nullptr),
    T_(M.T())
{
    if (M.v())
    {
        allocate();

        forAllBlockLinear(*this, i)
        {
            v_[i] = Zero;
        }
    }
}

template<class Type>
block<Type>::block(const block<Type>& M, const Type& v)
:
    refCount(),
    l_(M.l()),
    m_(M.m()),
    n_(M.n()),
    v_(nullptr),
    T_(M.T())
{
    if (M.v())
    {
        allocate();

        forAllBlockLinear(*this, i)
        {
            v_[i] = v;
        }
    }
}

template<class Type>
block<Type>::block(const tmp<block<Type>>& tM)
:
    refCount(),
    l_(tM->l()),
    m_(tM->m()),
    n_(tM->n()),
    v_(nullptr),
    T_(tM->T())
{
    if (tM.isTmp())
    {
        block<Type>& M = const_cast<block<Type>&>(tM());
        transfer(M);
    }
    else
    {
        allocate();
        *this = tM();
    }

    tM.clear();
}

template<class Type>
block<Type>::block(const tmp<block<Type>>& tM, const zero&)
:
    refCount(),
    l_(tM->l()),
    m_(tM->m()),
    n_(tM->n()),
    v_(nullptr),
    T_(tM->T())
{
    if (tM.isTmp())
    {
        block<Type>& M = const_cast<block<Type>&>(tM());
        transfer(M);
    }
    else
    {
        allocate();
    }

    *this = Zero;

    tM.clear();
}

template<class Type>
block<Type>::block(const tmp<block<Type>>& tM, const Type& v)
:
    refCount(),
    l_(tM->l()),
    m_(tM->m()),
    n_(tM->n()),
    v_(nullptr),
    T_(tM->T())
{
    if (tM.isTmp())
    {
        block<Type>& M = const_cast<block<Type>&>(tM());
        transfer(M);
    }
    else
    {
        allocate();
    }

    *this = v;

    tM.clear();
}

template<class Type>
block<Type>::block(const label reuse, block<Type>& M)
:
    refCount(),
    l_(M.l()),
    m_(M.m()),
    n_(M.n()),
    v_(nullptr),
    T_(M.T())
{
    if (reuse)
    {
        transfer(M);
    }
    else
    {
        allocate();
        *this = M;
    }
}

template<class Type>
block<Type>::block(const label reuse, block<Type>& M, const zero&)
:
    refCount(),
    l_(M.l()),
    m_(M.m()),
    n_(M.n()),
    v_(nullptr),
    T_(M.T())
{
    if (reuse)
    {
        transfer(M);
    }
    else
    {
        allocate();
    }

    *this = Zero;
}

template<class Type>
block<Type>::block(const label reuse, block<Type>& M, const Type& v)
:
    refCount(),
    l_(M.l()),
    m_(M.m()),
    n_(M.n()),
    v_(nullptr),
    T_(M.T())
{
    if (reuse)
    {
        transfer(M);
    }
    else
    {
        allocate();
    }

    *this = v;
}

template<class Type>
block<Type>::block(Istream& is)
:
    refCount(),
    l_(0),
    m_(0),
    n_(0),
    v_(nullptr),
    T_(eye)
{
    is >> *this;
}

// Destructor

template<class Type>
block<Type>::~block()
{
    if (v_)
    {
        delete[] v_;
    }
}

// Public

template<class Type>
void block<Type>::setSize
(
    const labelVector shape
)
{
    setSize(shape[0], shape[1], shape[2]);
}

template<class Type>
void block<Type>::setSize
(
    const int l,
    const int m,
    const int n
)
{
    block<Type> M(l, m, n, Zero);

    int minL = Foam::min(l, l_);
    int minM = Foam::min(m, m_);
    int minN = Foam::min(n, n_);

    for (int i = 0; i < minL; i++)
    {
        for (int j = 0; j < minM; j++)
        {
            for (int k = 0; k < minN; k++)
            {
                M(i,j,k) = (*this)(i,j,k);
            }
        }
    }

    transfer(M);
}

template<class Type>
void block<Type>::prepend
(
    const int d,
    const int s,
    const Type& v
)
{
    if (s <= 0)
    {
        FatalErrorInFunction
            << "Can only prepend a surface with positive thickness"
            << abort(FatalError);
    }

    if (d < 0 || d > 2)
    {
        FatalErrorInFunction
            << "Invalid prepend direction " << d
            << abort(FatalError);
    }

    const labelVector S(this->shape() + units[d]*s);

    block<Type> M(S,v);

    forAllBlock(*this, i, j, k)
    {
        M(labelVector(i,j,k)+units[d]*s) = this->operator()(i,j,k);
    }

    transfer(M);
}

template<class Type>
void block<Type>::append
(
    const int d,
    const int s,
    const Type& v
)
{
    if (s <= 0)
    {
        FatalErrorInFunction
            << "Can only append a surface with positive thickness"
            << abort(FatalError);
    }

    if (d < 0 || d > 2)
    {
        FatalErrorInFunction
            << "Invalid append direction " << d
            << abort(FatalError);
    }

    const labelVector S(this->shape() + units[d]*s);

    block<Type> M(S,v);

    forAllBlock(*this, i, j, k)
    {
        M(i,j,k) = this->operator()(i,j,k);
    }

    transfer(M);
}

template<class Type>
void block<Type>::transform(const labelTensor T)
{
    if (T == eye)
    {
        return;
    }

    if ((T.T() & T) != eye)
    {
        FatalErrorInFunction
            << "Tansformation matrix " << T << " is not orthogonal"
            << abort(FatalError);
    }

    // Transformed principle directions

    const labelVector x(T & unitX);
    const labelVector y(T & unitY);
    const labelVector z(T & unitZ);

    block<Type> M(cmptMag(T & this->shape()));

    // If the principle directions are negative, iteration should run from high
    // to low, requiring an offset

    const labelVector offset
    (
        ((x & unitXYZ) < 0 ? -x*(l_-1) : zeroXYZ)
      + ((y & unitXYZ) < 0 ? -y*(m_-1) : zeroXYZ)
      + ((z & unitXYZ) < 0 ? -z*(n_-1) : zeroXYZ)
    );

    // Iterate over all original points and construct the new index vectors from
    // the three transformed principle directions. The offset accounts for
    // negative principle directions. Some primitive types must be transformed
    // too. By default, pTransform does nothing.

    forAllBlock(*this, i, j, k)
    {
        const labelVector ijk(x*i + y*j + z*k + offset);

        M(ijk) = pTransform<Type>(T, this->operator()(i,j,k));
    }

    transfer(M);

    T_ = (T & T_);
}

template<class Type>
labelTensor block<Type>::reflect(const label a)
{
    if (a < 0 || a > 2)
    {
        FatalErrorInFunction
            << "Invalid axis when reflecting block" << endl
            << abort(FatalError);
    }

    const labelTensor T
    (
        a == 0 ? reflectX
      : a == 1 ? reflectY
      : reflectZ
    );

    transform(T);

    return T;
}

template<class Type>
labelTensor block<Type>::permute(const label a1, const label a2)
{
    if (a1 < 0 || a2 < 0 || a1 > 2 || a2 > 2)
    {
        FatalErrorInFunction
            << "Invalid axes specified for permutation" << endl
            << abort(FatalError);
    }

    if (a1 == a2)
    {
        return eye;
    }
    else
    {
        const labelTensor T
        (
            a1 == 0 && a2 == 1 ? permuteXY
          : a1 == 1 && a2 == 0 ? permuteXY
          : a1 == 0 && a2 == 2 ? permuteXZ
          : a1 == 2 && a2 == 0 ? permuteXZ
          : permuteYZ
        );

        transform(T);

        return T;
    }
}

template<class Type>
labelTensor block<Type>::rotate(const label st, const label a)
{
    if (a < 0 || a > 2)
    {
        FatalErrorInFunction
            << "Invalid axis specified for rotation" << endl
            << abort(FatalError);
    }

    const label s = (4+st%4)%4;

    if (s == 0)
    {
        return eye;
    }
    else
    {
        const labelTensor T
        (
            a == 0 ? (s == 1 ? rotateX1 : s == 2 ? rotateX2 : rotateX3)
          : a == 1 ? (s == 1 ? rotateY1 : s == 2 ? rotateY2 : rotateY3)
          :          (s == 1 ? rotateZ1 : s == 2 ? rotateZ2 : rotateZ3)
        );

        transform(T);

        return T;
    }
}

template<class Type>
void block<Type>::squeeze()
{
    if (l_ == 1 && m_ == 1 && n_ > 1)
    {
        l_ = n_;
        n_ = 1;
        T_ = (permuteXY&(permuteYZ&T_));
    }
    else if (l_ == 1 && m_ > 1 && n_ == 1)
    {
        l_ = m_;
        m_ = 1;
        T_ = (permuteXY&T_);
    }
    else if (l_ == 1 && m_ > 1 && n_ > 1)
    {
        l_ = m_;
        m_ = n_;
        n_ = 1;
        T_ = (permuteYZ&(permuteXY&T_));
    }
    else if (l_ > 1 && m_ == 1 && n_ > 1)
    {
        m_ = n_;
        n_ = 1;
        T_ = (permuteYZ&T_);
    }
}

template<class Type>
tmp<block<Type>> block<Type>::slicex(const label i, const bool sq) const
{
    const label ii = i < 0 ? l_+i : i;

    tmp<block<Type>> tM(new block<Type>(1,m_,n_));

    block<Type>& M = tM.ref();

    for (int j = 0; j < m_; j++)
    {
        for (int k = 0; k < n_; k++)
        {
            M(0,j,k) = this->operator()(ii,j,k);
        }
    }

    if (sq)
    {
        M.squeeze();
    }

    return tM;
}

template<class Type>
tmp<block<Type>> block<Type>::slicey(const label j, const bool sq) const
{
    const label jj = j < 0 ? m_+j : j;

    tmp<block<Type>> tM(new block<Type>(l_,1,n_));

    block<Type>& M = tM.ref();

    for (int i = 0; i < l_; i++)
    {
        for (int k = 0; k < n_; k++)
        {
            M(i,0,k) = this->operator()(i,jj,k);
        }
    }

    if (sq)
    {
        M.squeeze();
    }

    return tM;
}

template<class Type>
tmp<block<Type>> block<Type>::slicez(const label k, const bool sq) const
{
    const label kk = k < 0 ? n_+k : k;

    tmp<block<Type>> tM(new block<Type>(l_,m_,1));

    block<Type>& M = tM.ref();

    for (int i = 0; i < l_; i++)
    {
        for (int j = 0; j < m_; j++)
        {
            M(i,j,0) = this->operator()(i,j,kk);
        }
    }

    if (sq)
    {
        M.squeeze();
    }

    return tM;
}

template<class Type>
tmp<block<Type>> block<Type>::slice
(
    const label i,
    const label a,
    const bool sq
) const
{
    if (a < 0 || a > 2)
    {
        FatalErrorInFunction
            << "Invalid axis when slicing block" << endl
            << abort(FatalError);
    }

    if (a == 0)
    {
        return slicex(i, sq);
    }
    else if (a == 1)
    {
        return slicey(i, sq);
    }
    else
    {
        return slicez(i, sq);
    }
}

template<class Type>
tmp<block<Type>> block<Type>::sliceEdge(const label e, const bool sq) const
{
    const label side = (12+e%12)%12;
    const label dir = side/4;

    tmp<block<Type>> tM(new block<Type>(*this));

    block<Type>& M = tM.ref();

    if (dir == 0)
    {
        M.slice(-side%2,1);
        M.slice(-(side%4)/2,2);
    }
    if (dir == 1)
    {
        M.slice(-side%2,0);
        M.slice(-(side%4)/2,2);
    }
    else
    {
        M.slice(-side%2,0);
        M.slice(-(side%4)/2,1);
    }

    if (sq)
    {
        M.squeeze();
    }

    return tM;
}

template<class Type>
tmp<block<Type>> block<Type>::sliceVertex(const label v) const
{
    const label side = (8+v%8)%8;

    tmp<block<Type>> tM(new block<Type>(1,1,1));

    block<Type>& M = tM.ref();

    M(0,0,0) = (*this)(side%2,(side%4)/2,side/4);

    return tM;
}

template<class Type>
labelVector block<Type>::find(const Type& v) const
{
    forAllBlock(*this, i, j, k)
    {
        if ((*this)(i,j,k) == v)
        {
            return labelVector(i,j,k);
        }
    }

    return -unitXYZ;
}

template<class Type>
Type block<Type>::interp(const scalar i, const scalar j, const scalar k) const
{
    // Floor

    const label ii = i;
    const label jj = j;
    const label kk = k;

    // Set neighbor index as next higher value. For points exactly at a
    // face/edge/vertex, clip to the last point.

    const label iip = Foam::min(ii+1, l_-1);
    const label jjp = Foam::min(jj+1, m_-1);
    const label kkp = Foam::min(kk+1, n_-1);

    const scalar di = i-ii;
    const scalar dj = j-jj;
    const scalar dk = k-kk;

    return
        (1.0-di) * (1.0-dj) * (1.0-dk) * this->operator()(ii,  jj,  kk )
      + (    di) * (1.0-dj) * (1.0-dk) * this->operator()(iip, jj,  kk )
      + (1.0-di) * (    dj) * (1.0-dk) * this->operator()(ii,  jjp, kk )
      + (    di) * (    dj) * (1.0-dk) * this->operator()(iip, jjp, kk )
      + (1.0-di) * (1.0-dj) * (    dk) * this->operator()(ii,  jj,  kkp)
      + (    di) * (1.0-dj) * (    dk) * this->operator()(iip, jj,  kkp)
      + (1.0-di) * (    dj) * (    dk) * this->operator()(ii,  jjp, kkp)
      + (    di) * (    dj) * (    dk) * this->operator()(iip, jjp, kkp);
}

template<class Type>
Type block<Type>::interp(const vector ijk) const
{
    return this->interp(ijk.x(), ijk.y(), ijk.z());
}

template<class Type>
Istream& operator>>(Istream& is, block<Type>& M)
{
    List<List<List<Type>>> L;

    is >> L;

    const label l = L.size();
    const label m = L[0].size();
    const label n = L[0][0].size();

    forAll(L, i)
    {
        List<List<Type>>& Li = L[i];

        if (Li.size() != m)
        {
            FatalErrorInFunction
                << "Error while reading block from Istream: "
                << "provided block does not have uniform dimensions" << endl
                << abort(FatalError);
        }

        forAll(Li, j)
        {
            List<Type>& Lij = Li[j];

            if (Lij.size() != n)
            {
                FatalErrorInFunction
                    << "Error while reading block from Istream: "
                    << "provided block does not have uniform dimensions" << endl
                    << abort(FatalError);
            }
        }
    }

    M.clear();

    M.l_ = l;
    M.m_ = m;
    M.n_ = n;

    M.allocate();

    for (label i = 0; i < L.size(); i++)
    {
        List<List<Type>>& Li = L[i];

        for (label j = 0; j < Li.size(); j++)
        {
            List<Type>& Lij = L[i][j];

            for (label k = 0; k < Lij.size(); k++)
            {
                M(i,j,k) = Lij[k];
            }
        }
    }

    return is;
}

template<class Type>
Ostream& operator<<(Ostream& os, const block<Type>& M)
{
    List<List<List<Type>>> L(M.l_);

    forAll(L, i)
    {
        List<List<Type>>& Li = L[i];

        Li.setSize(M.m_);

        forAll(Li, j)
        {
            List<Type>& Lij = Li[j];

            Lij.setSize(M.n_);

            forAll(Lij, k)
            {
                Lij[k] = M(i,j,k);
            }
        }
    }

    os << L;

    return os;
}

}

}

#include "blockFunctions.C"
