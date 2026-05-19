namespace Foam
{

namespace briscola
{

template<class Type, int P>
inline void block<Type,P>::allocate()
{
    if (size_)
    {
        v_ = new Type[size_];

        #ifdef ZERO_INIT
        *this = Zero;
        #endif
    }
}

template<class Type, int P>
inline void block<Type,P>::clear()
{
    if (v_)
    {
        delete[] v_;
        v_ = nullptr;
    }

    l_ = 0;
    m_ = 0;
    n_ = 0;
    size_ = 0;
}

template<class Type, int P>
inline void block<Type,P>::reAllocate(int l, int m, int n)
{
    if (size_ != l*m*n)
    {
        clear();

        l_ = l;
        m_ = m;
        n_ = n;
        size_ = l*m*n;

        allocate();
    }
    else
    {
        // Same size but shape may be different

        l_ = l;
        m_ = m;
        n_ = n;
    }
}

template<class Type, int P>
inline void block<Type,P>::reAllocate(const labelVector& lmn)
{
    reAllocate(lmn.x(), lmn.y(), lmn.z());
}

template<class Type, int P>
inline void block<Type,P>::transfer(block<Type,P>& M)
{
    l_ = M.l_;
    m_ = M.m_;
    n_ = M.n_;
    size_ = l_*m_*n_;

    M.l_ = 0;
    M.m_ = 0;
    M.n_ = 0;
    M.size_ = 0;

    if (v_)
    {
        delete[] v_;
    }

    v_ = M.v_;

    M.v_ = nullptr;
}

template<class Type, int P>
inline void block<Type,P>::operator=(const block<Type,P>& M)
{
    if (this == &M)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    reAllocate(M.shape());

    if (v_)
        forAllBlockLinear(M, i)
            v_[i] = M.v()[i];
}

template<class Type, int P>
inline void block<Type,P>::operator=(const tmp<block<Type,P>>& tM)
{
    if (tM.isTmp())
    {
        block<Type,P>& M = const_cast<block<Type,P>&>(tM());
        transfer(M);
    }
    else
    {
        *this = tM();
    }

    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
inline void block<Type,P>::operator=(const Type& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] = v;
}

template<class Type, int P>
inline void block<Type,P>::operator=(const zero)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] = Zero;
}

template<class Type, int P>
inline void block<Type,P>::operator+=(const Type& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] += v;
}

template<class Type, int P>
inline void block<Type,P>::operator+=(const block<Type,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] += M.v()[i];
}

template<class Type, int P>
inline void block<Type,P>::operator+=(const tmp<block<Type,P>>& tM)
{
    *this += tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
inline void block<Type,P>::operator-=(const Type& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] -= v;
}

template<class Type, int P>
inline void block<Type,P>::operator-=(const block<Type,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] -= M.v()[i];
}

template<class Type, int P>
inline void block<Type,P>::operator-=(const tmp<block<Type,P>>& tM)
{
    *this -= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
inline void block<Type,P>::operator*=(const label& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] *= v;
}

template<class Type, int P>
inline void block<Type,P>::operator*=(const scalar& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] *= v;
}

template<class Type, int P>
inline void block<Type,P>::operator*=(const block<label,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] *= M.v()[i];
}

template<class Type, int P>
inline void block<Type,P>::operator*=(const tmp<block<label,P>>& tM)
{
    *this *= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
inline void block<Type,P>::operator*=(const block<scalar,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] *= M.v()[i];
}

template<class Type, int P>
inline void block<Type,P>::operator*=(const tmp<block<scalar,P>>& tM)
{
    *this *= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
inline void block<Type,P>::operator/=(const label& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] /= v;
}

template<class Type, int P>
inline void block<Type,P>::operator/=(const scalar& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] /= v;
}

template<class Type, int P>
inline void block<Type,P>::operator/=(const block<label,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] /= M.v()[i];
}

template<class Type, int P>
inline void block<Type,P>::operator/=(const tmp<block<label,P>>& tM)
{
    *this /= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
inline void block<Type,P>::operator/=(const block<scalar,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] /= M.v()[i];
}

template<class Type, int P>
inline void block<Type,P>::operator/=(const tmp<block<scalar,P>>& tM)
{
    *this /= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator=(const Type2& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] = v;
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator=(const block<Type2,P>& M)
{
    reAllocate(M.shape());

    if (v_)
        forAllBlockLinear(M, i)
            v_[i] = M.v()[i];
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator=(const tmp<block<Type2,P>>& tM)
{
    *this = tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator+=(const Type2& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] += v;
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator+=(const block<Type2,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] += M.v()[i];
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator+=(const tmp<block<Type2,P>>& tM)
{
    *this += tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator-=(const Type2& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] -= v;
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator-=(const block<Type2,P>& M)
{
    if (size_ != M.size())
    {
        FatalErrorInFunction
            << "Matrices do not have the same shape"
            << abort(FatalError);
    }

    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] -= M.v()[i];
}

template<class Type, int P>
template<class Type2>
inline void block<Type,P>::operator-=(const tmp<block<Type2,P>>& tM)
{
    *this -= tM();
    if (tM.isTmp())
        tM.clear();
}

}

}