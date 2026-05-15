namespace Foam
{

namespace briscola
{

template<class Type>
inline void block<Type>::allocate()
{
    if (size_)
    {
        v_ = new Type[size_];

        #ifdef ZERO_INIT
        *this = Zero;
        #endif
    }
}

template<class Type>
inline void block<Type>::clear()
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

template<class Type>
inline void block<Type>::reAllocate(int l, int m, int n)
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

template<class Type>
inline void block<Type>::reAllocate(const labelVector& lmn)
{
    reAllocate(lmn.x(), lmn.y(), lmn.z());
}

template<class Type>
inline void block<Type>::transfer(block<Type>& M)
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

template<class Type>
inline void block<Type>::operator=(const block<Type>& M)
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

template<class Type>
inline void block<Type>::operator=(const tmp<block<Type>>& tM)
{
    if (tM.isTmp())
    {
        block<Type>& M = const_cast<block<Type>&>(tM());
        transfer(M);
    }
    else
    {
        *this = tM();
    }

    if (tM.isTmp())
        tM.clear();
}

template<class Type>
inline void block<Type>::operator=(const Type& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] = v;
}

template<class Type>
inline void block<Type>::operator=(const zero)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] = Zero;
}

template<class Type>
inline void block<Type>::operator+=(const Type& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] += v;
}

template<class Type>
inline void block<Type>::operator+=(const block<Type>& M)
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

template<class Type>
inline void block<Type>::operator+=(const tmp<block<Type>>& tM)
{
    *this += tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
inline void block<Type>::operator-=(const Type& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] -= v;
}

template<class Type>
inline void block<Type>::operator-=(const block<Type>& M)
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

template<class Type>
inline void block<Type>::operator-=(const tmp<block<Type>>& tM)
{
    *this -= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
inline void block<Type>::operator*=(const label& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] *= v;
}

template<class Type>
inline void block<Type>::operator*=(const scalar& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] *= v;
}

template<class Type>
inline void block<Type>::operator*=(const block<label>& M)
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

template<class Type>
inline void block<Type>::operator*=(const tmp<block<label>>& tM)
{
    *this *= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
inline void block<Type>::operator*=(const block<scalar>& M)
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

template<class Type>
inline void block<Type>::operator*=(const tmp<block<scalar>>& tM)
{
    *this *= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
inline void block<Type>::operator/=(const label& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] /= v;
}

template<class Type>
inline void block<Type>::operator/=(const scalar& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] /= v;
}

template<class Type>
inline void block<Type>::operator/=(const block<label>& M)
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

template<class Type>
inline void block<Type>::operator/=(const tmp<block<label>>& tM)
{
    *this /= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
inline void block<Type>::operator/=(const block<scalar>& M)
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

template<class Type>
inline void block<Type>::operator/=(const tmp<block<scalar>>& tM)
{
    *this /= tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
template<class Type2>
inline void block<Type>::operator=(const Type2& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] = v;
}

template<class Type>
template<class Type2>
inline void block<Type>::operator=(const block<Type2>& M)
{
    reAllocate(M.shape());

    if (v_)
        forAllBlockLinear(M, i)
            v_[i] = M.v()[i];
}

template<class Type>
template<class Type2>
inline void block<Type>::operator=(const tmp<block<Type2>>& tM)
{
    *this = tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
template<class Type2>
inline void block<Type>::operator+=(const Type2& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] += v;
}

template<class Type>
template<class Type2>
inline void block<Type>::operator+=(const block<Type2>& M)
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

template<class Type>
template<class Type2>
inline void block<Type>::operator+=(const tmp<block<Type2>>& tM)
{
    *this += tM();
    if (tM.isTmp())
        tM.clear();
}

template<class Type>
template<class Type2>
inline void block<Type>::operator-=(const Type2& v)
{
    if (v_)
        forAllBlockLinear(*this, i)
            v_[i] -= v;
}

template<class Type>
template<class Type2>
inline void block<Type>::operator-=(const block<Type2>& M)
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

template<class Type>
template<class Type2>
inline void block<Type>::operator-=(const tmp<block<Type2>>& tM)
{
    *this -= tM();
    if (tM.isTmp())
        tM.clear();
}

}

}