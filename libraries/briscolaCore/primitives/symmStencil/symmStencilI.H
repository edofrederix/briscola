namespace Foam
{

inline symmStencil::symmStencil()
{}

inline symmStencil::symmStencil(const Foam::zero)
{
    this->v_[C] = Zero;
    this->v_[L] = Zero;
    this->v_[B] = Zero;
    this->v_[A] = Zero;
}

inline symmStencil::symmStencil
(
    const scalar& s1,
    const scalar& s2,
    const scalar& s3,
    const scalar& s4
)
{
    this->v_[C] = s1;
    this->v_[L] = s2;
    this->v_[B] = s3;
    this->v_[A] = s4;
}

inline symmStencil::symmStencil(const scalar& v)
{
    this->v_[C] = v;
    this->v_[L] = v;
    this->v_[B] = v;
    this->v_[A] = v;
}

inline symmStencil::symmStencil
(
    Istream& is
)
{
    is.readBegin("symmStencil");

    for (direction i=0; i<4; i++)
    {
        is >> v_[i];
    }

    is.readEnd("symmStencil");
    is.check("symmStencil::symmStencil(Istream&)");
}

inline direction symmStencil::size()
{
    return 4;
}

inline const scalar& symmStencil::component
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 4)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline scalar& symmStencil::component
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 4)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline void symmStencil::component
(
    scalar& s,
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 4)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    s = v_[d];
}

inline void symmStencil::replace
(
    const direction d,
    const scalar& s
)
{
    #ifdef FULLDEBUG
    if (d >= 4)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    v_[d] = s;
}

inline symmStencil symmStencil::uniform(const scalar& s)
{
    return symmStencil(s,s,s,s);
}

inline const scalar& Foam::symmStencil::center() const
{
    return this->v_[C];
}

inline const scalar& Foam::symmStencil::left() const
{
    return this->v_[L];
}

inline const scalar& Foam::symmStencil::bottom() const
{
    return this->v_[B];
}

inline const scalar& Foam::symmStencil::aft() const
{
    return this->v_[A];
}

inline scalar& Foam::symmStencil::center()
{
    return this->v_[C];
}

inline scalar& Foam::symmStencil::left()
{
    return this->v_[L];
}

inline scalar& Foam::symmStencil::bottom()
{
    return this->v_[B];
}

inline scalar& Foam::symmStencil::aft()
{
    return this->v_[A];
}

inline const scalar& symmStencil::operator[]
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 4)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline scalar& symmStencil::operator[]
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 4)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}

inline void symmStencil::operator+=
(
    const symmStencil& ss
)
{
    this->v_[C] += ss[C];
    this->v_[L] += ss[L];
    this->v_[B] += ss[B];
    this->v_[A] += ss[A];
}

inline void symmStencil::operator+=
(
    const diagStencil& st
)
{
    this->v_[C] += st[diagStencil::C];
}

inline void symmStencil::operator-=
(
    const symmStencil& ss
)
{
    this->v_[C] -= ss[C];
    this->v_[L] -= ss[L];
    this->v_[B] -= ss[B];
    this->v_[A] -= ss[A];
}

inline void symmStencil::operator-=
(
    const diagStencil& st
)
{
    this->v_[C] -= st[diagStencil::C];
}

inline void symmStencil::operator=
(
    const diagStencil& st
)
{
    this->v_[C] = st[diagStencil::C];
    this->v_[L] = Zero;
    this->v_[B] = Zero;
    this->v_[A] = Zero;
}

inline void symmStencil::operator=
(
    const scalar s
)
{
    this->v_[C] = s;
    this->v_[L] = s;
    this->v_[B] = s;
    this->v_[A] = s;
}

inline void symmStencil::operator=(const Foam::zero)
{
    this->v_[C] = Zero;
    this->v_[L] = Zero;
    this->v_[B] = Zero;
    this->v_[A] = Zero;
}

inline void symmStencil::operator*=
(
    const scalar s
)
{
    this->v_[C] *= s;
    this->v_[L] *= s;
    this->v_[B] *= s;
    this->v_[A] *= s;
}

inline void symmStencil::operator/=
(
    const scalar s
)
{
    this->v_[C] /= s;
    this->v_[L] /= s;
    this->v_[B] /= s;
    this->v_[A] /= s;
}

inline symmStencil operator+
(
    const symmStencil& st
)
{
    return st;
}


inline symmStencil operator-
(
    const symmStencil& st
)
{
    return symmStencil
    (
        -st[symmStencil::C],
        -st[symmStencil::L],
        -st[symmStencil::B],
        -st[symmStencil::A]
    );
}

inline symmStencil operator+
(
    const symmStencil& ss1,
    const symmStencil& ss2
)
{
    return symmStencil
    (
        ss1[symmStencil::C]+ss2[symmStencil::C],
        ss1[symmStencil::L]+ss2[symmStencil::L],
        ss1[symmStencil::B]+ss2[symmStencil::B],
        ss1[symmStencil::A]+ss2[symmStencil::A]
    );
}

inline symmStencil operator+
(
    const symmStencil& st1,
    const diagStencil& st2
)
{
    return symmStencil
    (
        st1[symmStencil::C]+st2[diagStencil::C],
        st1[symmStencil::L],
        st1[symmStencil::B],
        st1[symmStencil::A]
    );
}

inline symmStencil operator+
(
    const diagStencil& st1,
    const symmStencil& st2
)
{
    return st2+st1;
}

inline symmStencil operator-
(
    const symmStencil& ss1,
    const symmStencil& ss2
)
{
    return symmStencil
    (
        ss1[symmStencil::C]-ss2[symmStencil::C],
        ss1[symmStencil::L]-ss2[symmStencil::L],
        ss1[symmStencil::B]-ss2[symmStencil::B],
        ss1[symmStencil::A]-ss2[symmStencil::A]
    );
}

inline symmStencil operator-
(
    const symmStencil& st1,
    const diagStencil& st2
)
{
    return symmStencil
    (
        st1[symmStencil::C]-st2[diagStencil::C],
        st1[symmStencil::L],
        st1[symmStencil::B],
        st1[symmStencil::A]
    );
}

inline symmStencil operator-
(
    const diagStencil& st1,
    const symmStencil& st2
)
{
    return -(st2-st1);
}

inline symmStencil operator*
(
    scalar s,
    const symmStencil& ss
)
{
    return symmStencil
    (
        s*ss[symmStencil::C],
        s*ss[symmStencil::L],
        s*ss[symmStencil::B],
        s*ss[symmStencil::A]
    );
}

inline symmStencil operator*
(
    const symmStencil& ss,
    scalar s
)
{
    return s*ss;
}

inline symmStencil operator/
(
    const symmStencil& ss,
    scalar s
)
{
    return symmStencil
    (
        ss[symmStencil::C]/s,
        ss[symmStencil::L]/s,
        ss[symmStencil::B]/s,
        ss[symmStencil::A]/s
    );
}

inline bool operator==
(
    const symmStencil& st1,
    const symmStencil& st2
)
{
    return
        st1[symmStencil::C] == st2[symmStencil::C]
     && st1[symmStencil::L] == st2[symmStencil::L]
     && st1[symmStencil::B] == st2[symmStencil::B]
     && st1[symmStencil::A] == st2[symmStencil::A];
}

inline bool operator!=
(
    const symmStencil& st1,
    const symmStencil& st2
)
{
    return !(st1 == st2);
}

inline bool operator==
(
    const symmStencil& st1,
    const diagStencil& st2
)
{
    return
        st1[symmStencil::C] == st2[diagStencil::C]
     && st1[symmStencil::L] == 0.0
     && st1[symmStencil::B] == 0.0
     && st1[symmStencil::A] == 0.0;
}

inline bool operator==
(
    const diagStencil& st1,
    const symmStencil& st2
)
{
    return st2 == st1;
}

inline bool operator!=
(
    const symmStencil& st1,
    const diagStencil& st2
)
{
    return !(st1 == st2);
}

inline bool operator!=
(
    const diagStencil& st1,
    const symmStencil& st2
)
{
    return st2 != st1;
}

inline scalar stencilSum(const symmStencil& st)
{
    return
        st[symmStencil::C]
      + st[symmStencil::L]*2.0
      + st[symmStencil::B]*2.0
      + st[symmStencil::A]*2.0;
}

inline scalar neighborSum(const symmStencil& st)
{
    return
        st[symmStencil::L]*2.0
      + st[symmStencil::B]*2.0
      + st[symmStencil::A]*2.0;
}

inline scalar neighborSkewSum(const symmStencil& st)
{
    return 0.0;
}

template<>
class scalarProduct<scalar, symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class scalarProduct<symmStencil, scalar>
{
public:

    typedef symmStencil type;
};

template<>
class typeOfSum<symmStencil, diagStencil>
{
public:

    typedef symmStencil type;
};

template<>
class typeOfSum<diagStencil, symmStencil>
{
public:

    typedef symmStencil type;
};

}
