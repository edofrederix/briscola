#include "IOstreams.H"
#include <sstream>

namespace Foam
{

template<class Type>
inline Hex<Type>::Hex()
{}

template<class Type>
inline Hex<Type>::Hex(const Foam::zero)
{
    this->v_[L] = Zero;
    this->v_[R] = Zero;
    this->v_[B] = Zero;
    this->v_[T] = Zero;
    this->v_[A] = Zero;
    this->v_[F] = Zero;
}

template<class Type>
inline Hex<Type>::Hex
(
    const Type& v1,
    const Type& v2,
    const Type& v3,
    const Type& v4,
    const Type& v5,
    const Type& v6
)
{
    this->v_[L] = v1;
    this->v_[R] = v2;
    this->v_[B] = v3;
    this->v_[T] = v4;
    this->v_[A] = v5;
    this->v_[F] = v6;
}

template<class Type>
inline Hex<Type>::Hex
(
    Istream& is
)
{
    is.readBegin("Hex");

    for (direction i=0; i<6; i++)
    {
        is >> v_[i];
    }

    is.readEnd("Hex");
    is.check("Hex<Type>::Hex(Istream&)");
}

template<class Type>
inline const Type& Hex<Type>::component
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 6)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

template<class Type>
inline Type& Hex<Type>::component
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 6)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

template<class Type>
inline void Hex<Type>::component
(
    Type& c,
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 6)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    c = v_[d];
}

template<class Type>
inline void Hex<Type>::replace
(
    const direction d,
    const Type& c
)
{
    #ifdef FULLDEBUG
    if (d >= 6)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    v_[d] = c;
}

template<class Type>
inline Hex<Type> Hex<Type>::uniform(const Type& v)
{
    return Hex<Type>(v,v,v,v,v,v);
}

template<class Type>
inline const Type& Foam::Hex<Type>::left() const
{
    return this->v_[L];
}

template<class Type>
inline const Type& Foam::Hex<Type>::right() const
{
    return this->v_[R];
}

template<class Type>
inline const Type& Foam::Hex<Type>::bottom() const
{
    return this->v_[B];
}

template<class Type>
inline const Type& Foam::Hex<Type>::top() const
{
    return this->v_[T];
}

template<class Type>
inline const Type& Foam::Hex<Type>::aft() const
{
    return this->v_[A];
}

template<class Type>
inline const Type& Foam::Hex<Type>::fore() const
{
    return this->v_[F];
}

template<class Type>
inline Type& Foam::Hex<Type>::left()
{
    return this->v_[L];
}

template<class Type>
inline Type& Foam::Hex<Type>::right()
{
    return this->v_[R];
}

template<class Type>
inline Type& Foam::Hex<Type>::bottom()
{
    return this->v_[B];
}

template<class Type>
inline Type& Foam::Hex<Type>::top()
{
    return this->v_[T];
}

template<class Type>
inline Type& Foam::Hex<Type>::aft()
{
    return this->v_[A];
}

template<class Type>
inline Type& Foam::Hex<Type>::fore()
{
    return this->v_[F];
}

template<class Type>
inline Hex<Type> Hex<Type>::operator+() const
{
    return *this;
}

template<class Type>
inline Hex<Type> Hex<Type>::operator-() const
{
    return Hex<Type>
    (
      - this->v_[L],
      - this->v_[R],
      - this->v_[B],
      - this->v_[T],
      - this->v_[A],
      - this->v_[F]
    );
}

template<class Type>
inline const Type& Hex<Type>::operator[]
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 6)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

template<class Type>
inline Type& Hex<Type>::operator[]
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 6)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

template<class Type>
inline void Hex<Type>::operator+=
(
    const Hex<Type>& h
)
{
    this->v_[L] += h[L];
    this->v_[R] += h[R];
    this->v_[B] += h[B];
    this->v_[T] += h[T];
    this->v_[A] += h[A];
    this->v_[F] += h[F];
}

template<class Type>
inline void Hex<Type>::operator-=
(
    const Hex<Type>& h
)
{
    this->v_[L] -= h[L];
    this->v_[R] -= h[R];
    this->v_[B] -= h[B];
    this->v_[T] -= h[T];
    this->v_[A] -= h[A];
    this->v_[F] -= h[F];
}

template<class Type>
inline void Hex<Type>::operator=(const Foam::zero)
{
    this->v_[L] = Zero;
    this->v_[R] = Zero;
    this->v_[B] = Zero;
    this->v_[T] = Zero;
    this->v_[A] = Zero;
    this->v_[F] = Zero;
}

template<class Type>
inline void Hex<Type>::operator=(const Type& v)
{
    this->v_[L] = v;
    this->v_[R] = v;
    this->v_[B] = v;
    this->v_[T] = v;
    this->v_[A] = v;
    this->v_[F] = v;
}

template<class Type>
inline void Hex<Type>::operator*=
(
    const scalar s
)
{
    this->v_[L] *= s;
    this->v_[R] *= s;
    this->v_[B] *= s;
    this->v_[T] *= s;
    this->v_[A] *= s;
    this->v_[F] *= s;
}

template<class Type>
inline void Hex<Type>::operator*=
(
    const Hex<scalar>& h
)
{
    this->v_[L] *= h[L];
    this->v_[R] *= h[R];
    this->v_[B] *= h[B];
    this->v_[T] *= h[T];
    this->v_[A] *= h[A];
    this->v_[F] *= h[F];
}

template<class Type>
inline void Hex<Type>::operator/=
(
    const scalar s
)
{
    this->v_[L] /= s;
    this->v_[R] /= s;
    this->v_[B] /= s;
    this->v_[T] /= s;
    this->v_[A] /= s;
    this->v_[F] /= s;
}

template<class Type>
inline void Hex<Type>::operator/=
(
    const Hex<scalar>& h
)
{
    this->v_[L] /= h[L];
    this->v_[R] /= h[R];
    this->v_[B] /= h[B];
    this->v_[T] /= h[T];
    this->v_[A] /= h[A];
    this->v_[F] /= h[F];
}

// Define operator return types

template<class Type1, class Type2>
class typeOfSum<Type1,Hex<Type2>>
{
public:

    typedef Hex<typename typeOfSum<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class typeOfSum<Hex<Type1>,Type2>
{
public:

    typedef Hex<typename typeOfSum<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class typeOfSum<Hex<Type1>,Hex<Type2>>
{
public:

    typedef Hex<typename typeOfSum<Type1,Type2>::type> type;
};


template<class Type1, class Type2>
class outerProduct<Type1,Hex<Type2>>
{
public:

    typedef Hex<typename outerProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class outerProduct<Hex<Type1>,Type2>
{
public:

    typedef Hex<typename outerProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class outerProduct<Hex<Type1>,Hex<Type2>>
{
public:

    typedef Hex<typename outerProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class innerProduct<Type1,Hex<Type2>>
{
public:

    typedef Hex<typename innerProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class innerProduct<Hex<Type1>,Type2>
{
public:

    typedef Hex<typename innerProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class innerProduct<Hex<Type1>,Hex<Type2>>
{
public:

    typedef Hex<typename innerProduct<Type1,Type2>::type> type;
};


template<class Type1, class Type2>
class crossProduct<Type1,Hex<Type2>>
{
public:

    typedef Hex<typename crossProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class crossProduct<Hex<Type1>,Type2>
{
public:

    typedef Hex<typename crossProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class crossProduct<Hex<Type1>,Hex<Type2>>
{
public:

    typedef Hex<typename crossProduct<Type1,Type2>::type> type;
};


template<class Type1, class Type2>
class scalarProduct<Type1,Hex<Type2>>
{
public:

    typedef Hex<typename scalarProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class scalarProduct<Hex<Type1>,Type2>
{
public:

    typedef Hex<typename scalarProduct<Type1,Type2>::type> type;
};

template<class Type1, class Type2>
class scalarProduct<Hex<Type1>,Hex<Type2>>
{
public:

    typedef Hex<typename scalarProduct<Type1,Type2>::type> type;
};

// Summation/subtraction

template<class Type1, class Type2>
inline Hex<typename typeOfSum<Type1,Type2>::type> operator+
(
    const Hex<Type1>& h1,
    const Hex<Type2>& h2
)
{
    return Hex<typename typeOfSum<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] + h2[Hex<Type2>::L],
        h1[Hex<Type1>::R] + h2[Hex<Type2>::R],
        h1[Hex<Type1>::B] + h2[Hex<Type2>::B],
        h1[Hex<Type1>::T] + h2[Hex<Type2>::T],
        h1[Hex<Type1>::A] + h2[Hex<Type2>::A],
        h1[Hex<Type1>::F] + h2[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename typeOfSum<Type1,Type2>::type> operator-
(
    const Hex<Type1>& h1,
    const Hex<Type2>& h2
)
{
    return Hex<typename typeOfSum<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] - h2[Hex<Type2>::L],
        h1[Hex<Type1>::R] - h2[Hex<Type2>::R],
        h1[Hex<Type1>::B] - h2[Hex<Type2>::B],
        h1[Hex<Type1>::T] - h2[Hex<Type2>::T],
        h1[Hex<Type1>::A] - h2[Hex<Type2>::A],
        h1[Hex<Type1>::F] - h2[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename typeOfSum<Type1,Type2>::type> operator+
(
    const Hex<Type1>& h1,
    const Type2& v
)
{
    return Hex<typename typeOfSum<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] + v,
        h1[Hex<Type1>::R] + v,
        h1[Hex<Type1>::B] + v,
        h1[Hex<Type1>::T] + v,
        h1[Hex<Type1>::A] + v,
        h1[Hex<Type1>::F] + v
    );
}

template<class Type1, class Type2>
inline Hex<typename typeOfSum<Type1,Type2>::type> operator-
(
    const Hex<Type1>& h1,
    const Type2& v
)
{
    return Hex<typename typeOfSum<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] - v,
        h1[Hex<Type1>::R] - v,
        h1[Hex<Type1>::B] - v,
        h1[Hex<Type1>::T] - v,
        h1[Hex<Type1>::A] - v,
        h1[Hex<Type1>::F] - v
    );
}

template<class Type1, class Type2>
inline Hex<typename typeOfSum<Type1,Type2>::type> operator+
(
    const Type1& v,
    const Hex<Type2>& h2
)
{
    return Hex<typename typeOfSum<Type1,Type2>::type>
    (
        v + h2[Hex<Type2>::L],
        v + h2[Hex<Type2>::R],
        v + h2[Hex<Type2>::B],
        v + h2[Hex<Type2>::T],
        v + h2[Hex<Type2>::A],
        v + h2[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename typeOfSum<Type1,Type2>::type> operator-
(
    const Type1& v,
    const Hex<Type2>& h2
)
{
    return Hex<typename typeOfSum<Type1,Type2>::type>
    (
        v - h2[Hex<Type2>::L],
        v - h2[Hex<Type2>::R],
        v - h2[Hex<Type2>::B],
        v - h2[Hex<Type2>::T],
        v - h2[Hex<Type2>::A],
        v - h2[Hex<Type2>::F]
    );
}

// Multiplication

template<class Type>
inline Hex<Type> operator*
(
    const Hex<Type>& h1,
    const Hex<scalar>& h2
)
{
    return Hex<Type>
    (
        h1[Hex<Type>::L] * h2[Hex<scalar>::L],
        h1[Hex<Type>::R] * h2[Hex<scalar>::R],
        h1[Hex<Type>::B] * h2[Hex<scalar>::B],
        h1[Hex<Type>::T] * h2[Hex<scalar>::T],
        h1[Hex<Type>::A] * h2[Hex<scalar>::A],
        h1[Hex<Type>::F] * h2[Hex<scalar>::F]
    );
}

template<class Type>
inline Hex<Type> operator*
(
    const Hex<scalar>& h1,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        h1[Hex<scalar>::L] * h2[Hex<Type>::L],
        h1[Hex<scalar>::R] * h2[Hex<Type>::R],
        h1[Hex<scalar>::B] * h2[Hex<Type>::B],
        h1[Hex<scalar>::T] * h2[Hex<Type>::T],
        h1[Hex<scalar>::A] * h2[Hex<Type>::A],
        h1[Hex<scalar>::F] * h2[Hex<Type>::F]
    );
}

inline Hex<scalar> operator*
(
    const Hex<scalar>& h1,
    const Hex<scalar>& h2
)
{
    return Hex<scalar>
    (
        h1[Hex<scalar>::L] * h2[Hex<scalar>::L],
        h1[Hex<scalar>::R] * h2[Hex<scalar>::R],
        h1[Hex<scalar>::B] * h2[Hex<scalar>::B],
        h1[Hex<scalar>::T] * h2[Hex<scalar>::T],
        h1[Hex<scalar>::A] * h2[Hex<scalar>::A],
        h1[Hex<scalar>::F] * h2[Hex<scalar>::F]
    );
}

inline Hex<scalar> operator*
(
    const Hex<scalar>& h1,
    const scalar s
)
{
    return Hex<scalar>
    (
        h1[Hex<scalar>::L] * s,
        h1[Hex<scalar>::R] * s,
        h1[Hex<scalar>::B] * s,
        h1[Hex<scalar>::T] * s,
        h1[Hex<scalar>::A] * s,
        h1[Hex<scalar>::F] * s
    );
}

inline Hex<scalar> operator*
(
    const Hex<scalar>& h1,
    const label s
)
{
    return Hex<scalar>
    (
        h1[Hex<scalar>::L] * s,
        h1[Hex<scalar>::R] * s,
        h1[Hex<scalar>::B] * s,
        h1[Hex<scalar>::T] * s,
        h1[Hex<scalar>::A] * s,
        h1[Hex<scalar>::F] * s
    );
}

inline Hex<scalar> operator*
(
    const scalar s,
    const Hex<scalar>& h2
)
{
    return Hex<scalar>
    (
        s * h2[Hex<scalar>::L],
        s * h2[Hex<scalar>::R],
        s * h2[Hex<scalar>::B],
        s * h2[Hex<scalar>::T],
        s * h2[Hex<scalar>::A],
        s * h2[Hex<scalar>::F]
    );
}

inline Hex<scalar> operator*
(
    const label s,
    const Hex<scalar>& h2
)
{
    return Hex<scalar>
    (
        s * h2[Hex<scalar>::L],
        s * h2[Hex<scalar>::R],
        s * h2[Hex<scalar>::B],
        s * h2[Hex<scalar>::T],
        s * h2[Hex<scalar>::A],
        s * h2[Hex<scalar>::F]
    );
}

template<class Type>
inline Hex<Type> operator*
(
    const Hex<Type>& h1,
    const scalar s
)
{
    return Hex<Type>
    (
        h1[Hex<Type>::L] * s,
        h1[Hex<Type>::R] * s,
        h1[Hex<Type>::B] * s,
        h1[Hex<Type>::T] * s,
        h1[Hex<Type>::A] * s,
        h1[Hex<Type>::F] * s
    );
}

template<class Type>
inline Hex<Type> operator*
(
    const Hex<scalar>& h1,
    const Type& v
)
{
    return Hex<Type>
    (
        h1[Hex<Type>::L] * v,
        h1[Hex<Type>::R] * v,
        h1[Hex<Type>::B] * v,
        h1[Hex<Type>::T] * v,
        h1[Hex<Type>::A] * v,
        h1[Hex<Type>::F] * v
    );
}

template<class Type>
inline Hex<Type> operator*
(
    const scalar s,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        s * h2[Hex<Type>::L],
        s * h2[Hex<Type>::R],
        s * h2[Hex<Type>::B],
        s * h2[Hex<Type>::T],
        s * h2[Hex<Type>::A],
        s * h2[Hex<Type>::F]
    );
}

template<class Type>
inline Hex<Type> operator*
(
    const Type& v,
    const Hex<scalar>& h2
)
{
    return Hex<Type>
    (
        v * h2[Hex<scalar>::L],
        v * h2[Hex<scalar>::R],
        v * h2[Hex<scalar>::B],
        v * h2[Hex<scalar>::T],
        v * h2[Hex<scalar>::A],
        v * h2[Hex<scalar>::F]
    );
}

// Division

inline Hex<scalar> operator/
(
    const Hex<scalar>& h1,
    const Hex<scalar>& h2
)
{
    return Hex<scalar>
    (
        h1[Hex<scalar>::L] / h2[Hex<scalar>::L],
        h1[Hex<scalar>::R] / h2[Hex<scalar>::R],
        h1[Hex<scalar>::B] / h2[Hex<scalar>::B],
        h1[Hex<scalar>::T] / h2[Hex<scalar>::T],
        h1[Hex<scalar>::A] / h2[Hex<scalar>::A],
        h1[Hex<scalar>::F] / h2[Hex<scalar>::F]
    );
}

template<class Type>
inline Hex<Type> operator/
(
    const Hex<Type>& h1,
    const Hex<scalar>& h2
)
{
    return Hex<Type>
    (
        h1[Hex<Type>::L] / h2[Hex<scalar>::L],
        h1[Hex<Type>::R] / h2[Hex<scalar>::R],
        h1[Hex<Type>::B] / h2[Hex<scalar>::B],
        h1[Hex<Type>::T] / h2[Hex<scalar>::T],
        h1[Hex<Type>::A] / h2[Hex<scalar>::A],
        h1[Hex<Type>::F] / h2[Hex<scalar>::F]
    );
}

inline Hex<scalar> operator/
(
    const Hex<scalar>& h1,
    const scalar v
)
{
    return Hex<scalar>
    (
        h1[Hex<scalar>::L] / v,
        h1[Hex<scalar>::R] / v,
        h1[Hex<scalar>::B] / v,
        h1[Hex<scalar>::T] / v,
        h1[Hex<scalar>::A] / v,
        h1[Hex<scalar>::F] / v
    );
}

inline Hex<scalar> operator/
(
    const scalar v,
    const Hex<scalar>& h1
)
{
    return Hex<scalar>
    (
        v / h1[Hex<scalar>::L],
        v / h1[Hex<scalar>::R],
        v / h1[Hex<scalar>::B],
        v / h1[Hex<scalar>::T],
        v / h1[Hex<scalar>::A],
        v / h1[Hex<scalar>::F]
    );
}

template<class Type>
inline Hex<Type> operator/
(
    const Hex<Type>& h1,
    const scalar v
)
{
    return Hex<Type>
    (
        h1[Hex<Type>::L] / v,
        h1[Hex<Type>::R] / v,
        h1[Hex<Type>::B] / v,
        h1[Hex<Type>::T] / v,
        h1[Hex<Type>::A] / v,
        h1[Hex<Type>::F] / v
    );
}

template<class Type>
inline Hex<Type> operator/
(
    const Type& v,
    const Hex<scalar>& h2
)
{
    return Hex<Type>
    (
        v / h2[Hex<scalar>::L],
        v / h2[Hex<scalar>::R],
        v / h2[Hex<scalar>::B],
        v / h2[Hex<scalar>::T],
        v / h2[Hex<scalar>::A],
        v / h2[Hex<scalar>::F]
    );
}

// Outer products

template<class Type1, class Type2>
inline Hex<typename outerProduct<Type1,Type2>::type> operator*
(
    const Hex<Type1>& h1,
    const Hex<Type2>& h2
)
{
    return Hex<typename outerProduct<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] * h2[Hex<Type2>::L],
        h1[Hex<Type1>::R] * h2[Hex<Type2>::R],
        h1[Hex<Type1>::B] * h2[Hex<Type2>::B],
        h1[Hex<Type1>::T] * h2[Hex<Type2>::T],
        h1[Hex<Type1>::A] * h2[Hex<Type2>::A],
        h1[Hex<Type1>::F] * h2[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename outerProduct<Type1,Type2>::type> operator*
(
    const Type1& v,
    const Hex<Type2>& h
)
{
    return Hex<typename outerProduct<Type1,Type2>::type>
    (
        v * h[Hex<Type2>::L],
        v * h[Hex<Type2>::R],
        v * h[Hex<Type2>::B],
        v * h[Hex<Type2>::T],
        v * h[Hex<Type2>::A],
        v * h[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename outerProduct<Type1,Type2>::type> operator*
(
    const Hex<Type1>& h,
    const Type2& v
)
{
    return Hex<typename outerProduct<Type1,Type2>::type>
    (
        h[Hex<Type1>::L] * v,
        h[Hex<Type1>::R] * v,
        h[Hex<Type1>::B] * v,
        h[Hex<Type1>::T] * v,
        h[Hex<Type1>::A] * v,
        h[Hex<Type1>::F] * v
    );
}

// Scalar products

template<class Type1, class Type2>
inline Hex<typename scalarProduct<Type1,Type2>::type> operator&&
(
    const Hex<Type1>& h1,
    const Hex<Type2>& h2
)
{
    return Hex<typename scalarProduct<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] && h2[Hex<Type2>::L],
        h1[Hex<Type1>::R] && h2[Hex<Type2>::R],
        h1[Hex<Type1>::B] && h2[Hex<Type2>::B],
        h1[Hex<Type1>::T] && h2[Hex<Type2>::T],
        h1[Hex<Type1>::A] && h2[Hex<Type2>::A],
        h1[Hex<Type1>::F] && h2[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename scalarProduct<Type1,Type2>::type> operator&&
(
    const Hex<Type1>& h1,
    const Type2& v
)
{
    return Hex<typename scalarProduct<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] && v,
        h1[Hex<Type1>::R] && v,
        h1[Hex<Type1>::B] && v,
        h1[Hex<Type1>::T] && v,
        h1[Hex<Type1>::A] && v,
        h1[Hex<Type1>::F] && v
    );
}

template<class Type1, class Type2>
inline Hex<typename scalarProduct<Type1,Type2>::type> operator&&
(
    const Type1& v,
    const Hex<Type2>& h2
)
{
    return Hex<typename scalarProduct<Type1,Type2>::type>
    (
        v && h2[Hex<Type2>::L],
        v && h2[Hex<Type2>::R],
        v && h2[Hex<Type2>::B],
        v && h2[Hex<Type2>::T],
        v && h2[Hex<Type2>::A],
        v && h2[Hex<Type2>::F]
    );
}

// Inner products

template<class Type1, class Type2>
inline Hex<typename innerProduct<Type1,Type2>::type> operator&
(
    const Hex<Type1>& h1,
    const Hex<Type2>& h2
)
{
    return Hex<typename innerProduct<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] & h2[Hex<Type2>::L],
        h1[Hex<Type1>::R] & h2[Hex<Type2>::R],
        h1[Hex<Type1>::B] & h2[Hex<Type2>::B],
        h1[Hex<Type1>::T] & h2[Hex<Type2>::T],
        h1[Hex<Type1>::A] & h2[Hex<Type2>::A],
        h1[Hex<Type1>::F] & h2[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename innerProduct<Type1,Type2>::type> operator&
(
    const Hex<Type1>& h1,
    const Type2& v
)
{
    return Hex<typename innerProduct<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] & v,
        h1[Hex<Type1>::R] & v,
        h1[Hex<Type1>::B] & v,
        h1[Hex<Type1>::T] & v,
        h1[Hex<Type1>::A] & v,
        h1[Hex<Type1>::F] & v
    );
}

template<class Type1, class Type2>
inline Hex<typename innerProduct<Type1,Type2>::type> operator&
(
    const Type1& v,
    const Hex<Type2>& h2
)
{
    return Hex<typename innerProduct<Type1,Type2>::type>
    (
        v & h2[Hex<Type2>::L],
        v & h2[Hex<Type2>::R],
        v & h2[Hex<Type2>::B],
        v & h2[Hex<Type2>::T],
        v & h2[Hex<Type2>::A],
        v & h2[Hex<Type2>::F]
    );
}

// Cross products

template<class Type1, class Type2>
inline Hex<typename crossProduct<Type1,Type2>::type> operator^
(
    const Hex<Type1>& h1,
    const Hex<Type2>& h2
)
{
    return Hex<typename crossProduct<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] ^ h2[Hex<Type2>::L],
        h1[Hex<Type1>::R] ^ h2[Hex<Type2>::R],
        h1[Hex<Type1>::B] ^ h2[Hex<Type2>::B],
        h1[Hex<Type1>::T] ^ h2[Hex<Type2>::T],
        h1[Hex<Type1>::A] ^ h2[Hex<Type2>::A],
        h1[Hex<Type1>::F] ^ h2[Hex<Type2>::F]
    );
}

template<class Type1, class Type2>
inline Hex<typename crossProduct<Type1,Type2>::type> operator^
(
    const Hex<Type1>& h1,
    const Type2& v
)
{
    return Hex<typename crossProduct<Type1,Type2>::type>
    (
        h1[Hex<Type1>::L] ^ v,
        h1[Hex<Type1>::R] ^ v,
        h1[Hex<Type1>::B] ^ v,
        h1[Hex<Type1>::T] ^ v,
        h1[Hex<Type1>::A] ^ v,
        h1[Hex<Type1>::F] ^ v
    );
}

template<class Type1, class Type2>
inline Hex<typename crossProduct<Type1,Type2>::type> operator^
(
    const Type1& v,
    const Hex<Type2>& h2
)
{
    return Hex<typename crossProduct<Type1,Type2>::type>
    (
        v ^ h2[Hex<Type2>::L],
        v ^ h2[Hex<Type2>::R],
        v ^ h2[Hex<Type2>::B],
        v ^ h2[Hex<Type2>::T],
        v ^ h2[Hex<Type2>::A],
        v ^ h2[Hex<Type2>::F]
    );
}

// Comparison

template<class Type>
inline bool operator>
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    bool gt = true;
    for (direction i=0; i<6; ++i)
    {
        if (!(gt &= h1.v_[i] > h2.v_[i])) break;
    }
    return gt;
}

template<class Type>
inline bool operator<
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    bool lt = true;
    for (direction i=0; i<6; ++i)
    {
        if (!(lt &= h1.v_[i] < h2.v_[i])) break;
    }
    return lt;
}

template<class Type>
inline bool operator>=
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return !(h1 < h2);
}

template<class Type>
inline bool operator<=
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return !(h1 > h2);
}

template<class Type>
inline bool operator==
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return
        h1[Hex<Type>::L] == h2[Hex<Type>::L]
     && h1[Hex<Type>::R] == h2[Hex<Type>::R]
     && h1[Hex<Type>::B] == h2[Hex<Type>::B]
     && h1[Hex<Type>::T] == h2[Hex<Type>::T]
     && h1[Hex<Type>::A] == h2[Hex<Type>::A]
     && h1[Hex<Type>::F] == h2[Hex<Type>::F];
}

template<class Type>
inline bool operator!=
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return
        h1[Hex<Type>::L] != h2[Hex<Type>::L]
     || h1[Hex<Type>::R] != h2[Hex<Type>::R]
     || h1[Hex<Type>::B] != h2[Hex<Type>::B]
     || h1[Hex<Type>::T] != h2[Hex<Type>::T]
     || h1[Hex<Type>::A] != h2[Hex<Type>::A]
     || h1[Hex<Type>::F] != h2[Hex<Type>::F];
}

template<class Type>
inline Hex<scalar> magSqr
(
    const Hex<Type>& h
)
{
    return Hex<scalar>
    (
        magSqr(h[Hex<Type>::L]),
        magSqr(h[Hex<Type>::R]),
        magSqr(h[Hex<Type>::B]),
        magSqr(h[Hex<Type>::T]),
        magSqr(h[Hex<Type>::A]),
        magSqr(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<scalar> mag
(
    const Hex<Type>& h
)
{
    return Hex<scalar>
    (
        mag(h[Hex<Type>::L]),
        mag(h[Hex<Type>::R]),
        mag(h[Hex<Type>::B]),
        mag(h[Hex<Type>::T]),
        mag(h[Hex<Type>::A]),
        mag(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> cmptMultiply
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        cmptMultiply(h1[Hex<Type>::L], h2[Hex<Type>::L]),
        cmptMultiply(h1[Hex<Type>::R], h2[Hex<Type>::R]),
        cmptMultiply(h1[Hex<Type>::B], h2[Hex<Type>::B]),
        cmptMultiply(h1[Hex<Type>::T], h2[Hex<Type>::T]),
        cmptMultiply(h1[Hex<Type>::A], h2[Hex<Type>::A]),
        cmptMultiply(h1[Hex<Type>::F], h2[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> cmptPow
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        cmptPow(h1[Hex<Type>::L], h2[Hex<Type>::L]),
        cmptPow(h1[Hex<Type>::R], h2[Hex<Type>::R]),
        cmptPow(h1[Hex<Type>::B], h2[Hex<Type>::B]),
        cmptPow(h1[Hex<Type>::T], h2[Hex<Type>::T]),
        cmptPow(h1[Hex<Type>::A], h2[Hex<Type>::A]),
        cmptPow(h1[Hex<Type>::F], h2[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> cmptDivide
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        cmptDivide(h1[Hex<Type>::L], h2[Hex<Type>::L]),
        cmptDivide(h1[Hex<Type>::R], h2[Hex<Type>::R]),
        cmptDivide(h1[Hex<Type>::B], h2[Hex<Type>::B]),
        cmptDivide(h1[Hex<Type>::T], h2[Hex<Type>::T]),
        cmptDivide(h1[Hex<Type>::A], h2[Hex<Type>::A]),
        cmptDivide(h1[Hex<Type>::F], h2[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<typename pTraits<Type>::cmptType> cmptMax
(
    const Hex<Type>& h
)
{
    return Hex<typename pTraits<Type>::cmptType>
    (
        cmptMax(h[Hex<Type>::L]),
        cmptMax(h[Hex<Type>::R]),
        cmptMax(h[Hex<Type>::B]),
        cmptMax(h[Hex<Type>::T]),
        cmptMax(h[Hex<Type>::A]),
        cmptMax(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<typename pTraits<Type>::cmptType> cmptMin
(
    const Hex<Type>& h
)
{
    return Hex<typename pTraits<Type>::cmptType>
    (
        cmptMin(h[Hex<Type>::L]),
        cmptMin(h[Hex<Type>::R]),
        cmptMin(h[Hex<Type>::B]),
        cmptMin(h[Hex<Type>::T]),
        cmptMin(h[Hex<Type>::A]),
        cmptMin(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<typename pTraits<Type>::cmptType> cmptSum
(
    const Hex<Type>& h
)
{
    return Hex<typename pTraits<Type>::cmptType>
    (
        cmptSum(h[Hex<Type>::L]),
        cmptSum(h[Hex<Type>::R]),
        cmptSum(h[Hex<Type>::B]),
        cmptSum(h[Hex<Type>::T]),
        cmptSum(h[Hex<Type>::A]),
        cmptSum(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<typename pTraits<Type>::cmptType> cmptAv
(
    const Hex<Type>& h
)
{
    return cmptSum(h)/pTraits<Type>::nComponents;
}

template<class Type>
inline Hex<typename pTraits<Type>::cmptType> cmptProduct
(
    const Hex<Type>& h
)
{
    return Hex<typename pTraits<Type>::cmptType>
    (
        cmptProduct(h[Hex<Type>::L]),
        cmptProduct(h[Hex<Type>::R]),
        cmptProduct(h[Hex<Type>::B]),
        cmptProduct(h[Hex<Type>::T]),
        cmptProduct(h[Hex<Type>::A]),
        cmptProduct(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> cmptSqr
(
    const Hex<Type>& h
)
{
    return Hex<Type>
    (
        cmptSqr(h[Hex<Type>::L]),
        cmptSqr(h[Hex<Type>::R]),
        cmptSqr(h[Hex<Type>::B]),
        cmptSqr(h[Hex<Type>::T]),
        cmptSqr(h[Hex<Type>::A]),
        cmptSqr(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> cmptMag
(
    const Hex<Type>& h
)
{
    return Hex<Type>
    (
        cmptMag(h[Hex<Type>::L]),
        cmptMag(h[Hex<Type>::R]),
        cmptMag(h[Hex<Type>::B]),
        cmptMag(h[Hex<Type>::T]),
        cmptMag(h[Hex<Type>::A]),
        cmptMag(h[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> max
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        max(h1[Hex<Type>::L], h2[Hex<Type>::L]),
        max(h1[Hex<Type>::R], h2[Hex<Type>::R]),
        max(h1[Hex<Type>::B], h2[Hex<Type>::B]),
        max(h1[Hex<Type>::T], h2[Hex<Type>::T]),
        max(h1[Hex<Type>::A], h2[Hex<Type>::A]),
        max(h1[Hex<Type>::F], h2[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> min
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        max(h1[Hex<Type>::L], h2[Hex<Type>::L]),
        max(h1[Hex<Type>::R], h2[Hex<Type>::R]),
        max(h1[Hex<Type>::B], h2[Hex<Type>::B]),
        max(h1[Hex<Type>::T], h2[Hex<Type>::T]),
        max(h1[Hex<Type>::A], h2[Hex<Type>::A]),
        max(h1[Hex<Type>::F], h2[Hex<Type>::F])
    );
}

template<class Type>
inline Hex<Type> minMod
(
    const Hex<Type>& h1,
    const Hex<Type>& h2
)
{
    return Hex<Type>
    (
        minMod(h1[Hex<Type>::L], h2[Hex<Type>::L]),
        minMod(h1[Hex<Type>::R], h2[Hex<Type>::R]),
        minMod(h1[Hex<Type>::B], h2[Hex<Type>::B]),
        minMod(h1[Hex<Type>::T], h2[Hex<Type>::T]),
        minMod(h1[Hex<Type>::A], h2[Hex<Type>::A]),
        minMod(h1[Hex<Type>::F], h2[Hex<Type>::F])
    );
}

template<class Type>
inline Type neighborSum
(
    const Hex<Type>& h
)
{
    return
        h[Hex<Type>::L]
      + h[Hex<Type>::R]
      + h[Hex<Type>::B]
      + h[Hex<Type>::T]
      + h[Hex<Type>::A]
      + h[Hex<Type>::F];
}

template<class Type>
inline Type neighborSkewSum
(
    const Hex<Type>& h
)
{
    return
      - h[Hex<Type>::L]
      + h[Hex<Type>::R]
      - h[Hex<Type>::B]
      + h[Hex<Type>::T]
      - h[Hex<Type>::A]
      + h[Hex<Type>::F];
}

template<class Type>
inline Hex<Type> dot
(
    const scalar s,
    const Hex<Type>& h2
)
{
    return s*h2;
}

template<class Type>
inline Hex<Type> dot
(
    const Hex<Type>& h1,
    const scalar s
)
{
    return h1*s;
}

template<class Type>
inline Hex<Type> dot
(
    const Hex<scalar>& h1,
    const Hex<Type>& h2
)
{
    return h1*h2;
}

template<class Type>
inline Hex<Type> dot
(
    const Hex<Type>& h1,
    const Hex<scalar>& h2
)
{
    return h1*h2;
}

template<class Type1, class Type2>
inline Hex<typename innerProduct<Type1,Type2>::type> dot
(
    const Hex<Type1>& h1,
    const Hex<Type2>& h2
)
{
    return h1 & h2;
}

template<class Type1, class Type2>
inline Hex<typename innerProduct<Type1,Type2>::type> dot
(
    const Type1& v,
    const Hex<Type2>& h2
)
{
    return v & h2;
}

template<class Type1, class Type2>
inline Hex<typename innerProduct<Type1,Type2>::type> dot
(
    const Hex<Type1>& h1,
    const Type2& v
)
{
    return h1 & v;
}

template<class Type>
inline direction findMax
(
    const Hex<Type>& h
)
{
    direction index = 0;
    for (direction i=1; i<6; ++i)
    {
        index = h[index] > h[i] ? index : i;
    }
    return index;
}

template<class Type>
inline direction findMin
(
    const Hex<Type>& h
)
{
    direction index = 0;
    for (direction i=1; i<6; ++i)
    {
        index = h[index] < h[i] ? index : i;
    }
    return index;
}


}
