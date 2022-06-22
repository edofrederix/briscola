namespace Foam
{

inline stencil::stencil()
{}

inline stencil::stencil(const Foam::zero)
{
    this->v_[C] = Zero;
    this->v_[L] = Zero;
    this->v_[R] = Zero;
    this->v_[B] = Zero;
    this->v_[T] = Zero;
    this->v_[A] = Zero;
    this->v_[F] = Zero;
}

inline stencil::stencil
(
    const scalar& s1,
    const scalar& s2,
    const scalar& s3,
    const scalar& s4,
    const scalar& s5,
    const scalar& s6,
    const scalar& s7
)
{
    this->v_[C] = s1;
    this->v_[L] = s2;
    this->v_[R] = s3;
    this->v_[B] = s4;
    this->v_[T] = s5;
    this->v_[A] = s6;
    this->v_[F] = s7;
}

inline stencil::stencil(const scalar& v)
{
    this->v_[C] = v;
    this->v_[L] = v;
    this->v_[R] = v;
    this->v_[B] = v;
    this->v_[T] = v;
    this->v_[A] = v;
    this->v_[F] = v;
}

inline stencil::stencil(const stencil& st)
{
    this->v_[C] = st[C];
    this->v_[L] = st[L];
    this->v_[R] = st[R];
    this->v_[B] = st[B];
    this->v_[T] = st[T];
    this->v_[A] = st[A];
    this->v_[F] = st[F];
}

inline stencil::stencil(const symmStencil& st)
{
    this->v_[C] = st[symmStencil::C];
    this->v_[L] = st[symmStencil::L];
    this->v_[R] = st[symmStencil::L];
    this->v_[B] = st[symmStencil::B];
    this->v_[T] = st[symmStencil::B];
    this->v_[A] = st[symmStencil::A];
    this->v_[F] = st[symmStencil::A];
}

inline stencil::stencil(const diagStencil& st)
{
    this->v_[C] = st[diagStencil::C];
    this->v_[L] = Zero;
    this->v_[R] = Zero;
    this->v_[B] = Zero;
    this->v_[T] = Zero;
    this->v_[A] = Zero;
    this->v_[F] = Zero;
}

inline stencil::stencil
(
    Istream& is
)
{
    is.readBegin("stencil");

    for (direction i=0; i<7; i++)
    {
        is >> v_[i];
    }

    is.readEnd("stencil");
    is.check("stencil::stencil(Istream&)");
}

inline direction stencil::size()
{
    return 7;
}

inline const scalar& stencil::component
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 7)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline scalar& stencil::component
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 7)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline void stencil::component
(
    scalar& s,
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 7)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    s = v_[d];
}

inline void stencil::replace
(
    const direction d,
    const scalar& s
)
{
    #ifdef FULLDEBUG
    if (d >= 7)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    v_[d] = s;
}

inline stencil stencil::uniform(const scalar& s)
{
    return stencil(s,s,s,s,s,s,s);
}

inline const scalar& Foam::stencil::center() const
{
    return this->v_[C];
}

inline const scalar& Foam::stencil::left() const
{
    return this->v_[L];
}

inline const scalar& Foam::stencil::right() const
{
    return this->v_[R];
}

inline const scalar& Foam::stencil::bottom() const
{
    return this->v_[B];
}

inline const scalar& Foam::stencil::top() const
{
    return this->v_[T];
}

inline const scalar& Foam::stencil::aft() const
{
    return this->v_[A];
}

inline const scalar& Foam::stencil::fore() const
{
    return this->v_[F];
}

inline scalar& Foam::stencil::center()
{
    return this->v_[C];
}

inline scalar& Foam::stencil::left()
{
    return this->v_[L];
}

inline scalar& Foam::stencil::right()
{
    return this->v_[R];
}

inline scalar& Foam::stencil::bottom()
{
    return this->v_[B];
}

inline scalar& Foam::stencil::top()
{
    return this->v_[T];
}

inline scalar& Foam::stencil::aft()
{
    return this->v_[A];
}

inline scalar& Foam::stencil::fore()
{
    return this->v_[F];
}

inline const scalar& stencil::operator[]
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 7)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline scalar& stencil::operator[]
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 7)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline void stencil::operator+=
(
    const stencil& st
)
{
    this->v_[C] += st[C];
    this->v_[L] += st[L];
    this->v_[R] += st[R];
    this->v_[B] += st[B];
    this->v_[T] += st[T];
    this->v_[A] += st[A];
    this->v_[F] += st[F];
}

inline void stencil::operator-=
(
    const stencil& st
)
{
    this->v_[C] -= st[C];
    this->v_[L] -= st[L];
    this->v_[R] -= st[R];
    this->v_[B] -= st[B];
    this->v_[T] -= st[T];
    this->v_[A] -= st[A];
    this->v_[F] -= st[F];
}

inline void stencil::operator+=
(
    const symmStencil& st
)
{
    this->v_[C] += st[symmStencil::C];
    this->v_[L] += st[symmStencil::L];
    this->v_[R] += st[symmStencil::L];
    this->v_[B] += st[symmStencil::B];
    this->v_[T] += st[symmStencil::B];
    this->v_[A] += st[symmStencil::A];
    this->v_[F] += st[symmStencil::A];
}

inline void stencil::operator+=
(
    const diagStencil& st
)
{
    this->v_[C] += st[diagStencil::C];
}

inline void stencil::operator-=
(
    const symmStencil& st
)
{
    this->v_[C] -= st[symmStencil::C];
    this->v_[L] -= st[symmStencil::L];
    this->v_[R] -= st[symmStencil::L];
    this->v_[B] -= st[symmStencil::B];
    this->v_[T] -= st[symmStencil::B];
    this->v_[A] -= st[symmStencil::A];
    this->v_[F] -= st[symmStencil::A];
}

inline void stencil::operator-=
(
    const diagStencil& st
)
{
    this->v_[C] -= st[diagStencil::C];
}

inline void stencil::operator=(const Foam::zero)
{
    this->v_[C] = Zero;
    this->v_[L] = Zero;
    this->v_[R] = Zero;
    this->v_[B] = Zero;
    this->v_[T] = Zero;
    this->v_[A] = Zero;
    this->v_[F] = Zero;
}

inline void stencil::operator=
(
    const scalar s
)
{
    this->v_[C] = s;
    this->v_[L] = s;
    this->v_[R] = s;
    this->v_[B] = s;
    this->v_[T] = s;
    this->v_[A] = s;
    this->v_[F] = s;
}

inline void stencil::operator=
(
    const stencil& st
)
{
    this->v_[C] = st[C];
    this->v_[L] = st[L];
    this->v_[R] = st[R];
    this->v_[B] = st[B];
    this->v_[T] = st[T];
    this->v_[A] = st[A];
    this->v_[F] = st[F];
}

inline void stencil::operator=
(
    const symmStencil& st
)
{
    this->v_[C] = st[symmStencil::C];
    this->v_[L] = st[symmStencil::L];
    this->v_[R] = st[symmStencil::L];
    this->v_[B] = st[symmStencil::B];
    this->v_[T] = st[symmStencil::B];
    this->v_[A] = st[symmStencil::A];
    this->v_[F] = st[symmStencil::A];
}

inline void stencil::operator=
(
    const diagStencil& st
)
{
    this->v_[C] = st[diagStencil::C];
    this->v_[L] = Zero;
    this->v_[R] = Zero;
    this->v_[B] = Zero;
    this->v_[T] = Zero;
    this->v_[A] = Zero;
    this->v_[F] = Zero;
}

inline void stencil::operator=
(
    const hexScalar& hs
)
{
    this->v_[C] = Zero;
    this->v_[L] = hs[hexScalar::L];
    this->v_[R] = hs[hexScalar::R];
    this->v_[B] = hs[hexScalar::B];
    this->v_[T] = hs[hexScalar::T];
    this->v_[A] = hs[hexScalar::A];
    this->v_[F] = hs[hexScalar::F];
}

inline void stencil::operator*=
(
    const scalar s
)
{
    this->v_[C] *= s;
    this->v_[L] *= s;
    this->v_[R] *= s;
    this->v_[B] *= s;
    this->v_[T] *= s;
    this->v_[A] *= s;
    this->v_[F] *= s;
}

inline void stencil::operator/=
(
    const scalar s
)
{
    this->v_[C] /= s;
    this->v_[L] /= s;
    this->v_[R] /= s;
    this->v_[B] /= s;
    this->v_[T] /= s;
    this->v_[A] /= s;
    this->v_[F] /= s;
}

inline stencil operator+
(
    const stencil& st
)
{
    return st;
}

inline stencil operator-
(
    const stencil& st
)
{
    return stencil
    (
        -st[stencil::C],
        -st[stencil::L],
        -st[stencil::R],
        -st[stencil::B],
        -st[stencil::T],
        -st[stencil::A],
        -st[stencil::F]
    );
}

inline stencil operator+
(
    const stencil& st1,
    const stencil& st2
)
{
    return stencil
    (
        st1[stencil::C]+st2[stencil::C],
        st1[stencil::L]+st2[stencil::L],
        st1[stencil::R]+st2[stencil::R],
        st1[stencil::B]+st2[stencil::B],
        st1[stencil::T]+st2[stencil::T],
        st1[stencil::A]+st2[stencil::A],
        st1[stencil::F]+st2[stencil::F]
    );
}

inline stencil operator+
(
    const stencil& st1,
    const symmStencil& st2
)
{
    return stencil
    (
        st1[stencil::C]+st2[symmStencil::C],
        st1[stencil::L]+st2[symmStencil::L],
        st1[stencil::R]+st2[symmStencil::L],
        st1[stencil::B]+st2[symmStencil::B],
        st1[stencil::T]+st2[symmStencil::B],
        st1[stencil::A]+st2[symmStencil::A],
        st1[stencil::F]+st2[symmStencil::A]
    );
}

inline stencil operator+
(
    const stencil& st1,
    const diagStencil& st2
)
{
    return stencil
    (
        st1[stencil::C]+st2[diagStencil::C],
        st1[stencil::L],
        st1[stencil::R],
        st1[stencil::B],
        st1[stencil::T],
        st1[stencil::A],
        st1[stencil::F]
    );
}

inline stencil operator+
(
    const symmStencil& st1,
    const stencil& st2
)
{
    return st2+st1;
}

inline stencil operator+
(
    const diagStencil& st1,
    const stencil& st2
)
{
    return st2+st1;
}

inline stencil operator-
(
    const stencil& st1,
    const stencil& st2
)
{
    return stencil
    (
        st1[stencil::C]-st2[stencil::C],
        st1[stencil::L]-st2[stencil::L],
        st1[stencil::R]-st2[stencil::R],
        st1[stencil::B]-st2[stencil::B],
        st1[stencil::T]-st2[stencil::T],
        st1[stencil::A]-st2[stencil::A],
        st1[stencil::F]-st2[stencil::F]
    );
}

inline stencil operator-
(
    const stencil& st1,
    const symmStencil& st2
)
{
    return stencil
    (
        st1[stencil::C]-st2[symmStencil::C],
        st1[stencil::L]-st2[symmStencil::L],
        st1[stencil::R]-st2[symmStencil::L],
        st1[stencil::B]-st2[symmStencil::B],
        st1[stencil::T]-st2[symmStencil::B],
        st1[stencil::A]-st2[symmStencil::A],
        st1[stencil::F]-st2[symmStencil::A]
    );
}

inline stencil operator-
(
    const stencil& st1,
    const diagStencil& st2
)
{
    return stencil
    (
        st1[stencil::C]-st2[diagStencil::C],
        st1[stencil::L],
        st1[stencil::R],
        st1[stencil::B],
        st1[stencil::T],
        st1[stencil::A],
        st1[stencil::F]
    );
}

inline stencil operator-
(
    const symmStencil& st1,
    const stencil& st2
)
{
    return -(st2-st1);
}

inline stencil operator-
(
    const diagStencil& st1,
    const stencil& st2
)
{
    return -(st2-st1);
}

inline stencil operator*
(
    scalar s,
    const stencil& st
)
{
    return stencil
    (
        s*st[stencil::C],
        s*st[stencil::L],
        s*st[stencil::R],
        s*st[stencil::B],
        s*st[stencil::T],
        s*st[stencil::A],
        s*st[stencil::F]
    );
}

inline stencil operator*
(
    const stencil& st,
    scalar s
)
{
    return s*st;
}

inline stencil operator/
(
    const stencil& st,
    scalar s
)
{
    return stencil
    (
        st[stencil::C]/s,
        st[stencil::L]/s,
        st[stencil::R]/s,
        st[stencil::B]/s,
        st[stencil::T]/s,
        st[stencil::A]/s,
        st[stencil::F]/s
    );
}

inline bool operator==
(
    const stencil& st1,
    const stencil& st2
)
{
    return
        st1[stencil::C] == st2[stencil::C]
     && st1[stencil::L] == st2[stencil::L]
     && st1[stencil::R] == st2[stencil::R]
     && st1[stencil::B] == st2[stencil::B]
     && st1[stencil::T] == st2[stencil::T]
     && st1[stencil::A] == st2[stencil::A]
     && st1[stencil::F] == st2[stencil::F];
}

inline bool operator!=
(
    const stencil& st1,
    const stencil& st2
)
{
    return !(st1 == st2);
}

inline bool operator==
(
    const stencil& st1,
    const symmStencil& st2
)
{
    return
        st1[stencil::C] == st2[symmStencil::C]
     && st1[stencil::L] == st2[symmStencil::L]
     && st1[stencil::R] == st2[symmStencil::L]
     && st1[stencil::B] == st2[symmStencil::B]
     && st1[stencil::T] == st2[symmStencil::B]
     && st1[stencil::A] == st2[symmStencil::A]
     && st1[stencil::F] == st2[symmStencil::A];
}

inline bool operator==
(
    const symmStencil& st1,
    const stencil& st2
)
{
    return st2 == st1;
}

inline bool operator!=
(
    const stencil& st1,
    const symmStencil& st2
)
{
    return !(st1 == st2);
}

inline bool operator!=
(
    const symmStencil& st1,
    const stencil& st2
)
{
    return st2 != st1;
}

inline bool operator==
(
    const stencil& st1,
    const diagStencil& st2
)
{
    return
        st1[stencil::C] == st2[diagStencil::C]
     && st1[stencil::L] == 0.0
     && st1[stencil::R] == 0.0
     && st1[stencil::B] == 0.0
     && st1[stencil::T] == 0.0
     && st1[stencil::A] == 0.0
     && st1[stencil::F] == 0.0;
}

inline bool operator==
(
    const diagStencil& st1,
    const stencil& st2
)
{
    return st2 == st1;
}

inline bool operator!=
(
    const stencil& st1,
    const diagStencil& st2
)
{
    return !(st1 == st2);
}

inline bool operator!=
(
    const diagStencil& st1,
    const stencil& st2
)
{
    return st2 != st1;
}

inline scalar stencilSum(const stencil& st)
{
    return
        st[stencil::C]
      + st[stencil::L]
      + st[stencil::R]
      + st[stencil::B]
      + st[stencil::T]
      + st[stencil::A]
      + st[stencil::F];
}

inline scalar neighborSum(const stencil& st)
{
    return
        st[stencil::L]
      + st[stencil::R]
      + st[stencil::B]
      + st[stencil::T]
      + st[stencil::A]
      + st[stencil::F];
}

inline scalar neighborSkewSum(const stencil& st)
{
    return
      - st[stencil::L]
      + st[stencil::R]
      - st[stencil::B]
      + st[stencil::T]
      - st[stencil::A]
      + st[stencil::F];
}

template<>
class scalarProduct<scalar, stencil>
{
public:

    typedef stencil type;
};

template<>
class scalarProduct<stencil, scalar>
{
public:

    typedef stencil type;
};

template<>
class typeOfSum<stencil, symmStencil>
{
public:

    typedef stencil type;
};

template<>
class typeOfSum<symmStencil, stencil>
{
public:

    typedef stencil type;
};

template<>
class typeOfSum<stencil, diagStencil>
{
public:

    typedef stencil type;
};

template<>
class typeOfSum<diagStencil, stencil>
{
public:

    typedef stencil type;
};

}
