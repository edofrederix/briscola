namespace Foam
{

inline diagStencil::diagStencil()
{}

inline diagStencil::diagStencil(const Foam::zero)
{
    this->v_ = Zero;
}

inline diagStencil::diagStencil
(
    const scalar& s1
)
{
    this->v_ = s1;
}

inline diagStencil::diagStencil
(
    Istream& is
)
{
    is.readBegin("diagStencil");

    is >> v_;

    is.readEnd("diagStencil");
    is.check("diagStencil::diagStencil(Istream&)");
}

inline direction diagStencil::size()
{
    return 1;
}

inline const scalar& diagStencil::component
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 1)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_;
}

inline scalar& diagStencil::component
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 1)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_;
}

inline void diagStencil::component
(
    scalar& s,
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 1)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    s = v_;
}

inline void diagStencil::replace
(
    const direction d,
    const scalar& s
)
{
    #ifdef FULLDEBUG
    if (d >= 1)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    v_ = s;
}

inline diagStencil diagStencil::uniform(const scalar& s)
{
    return diagStencil(s);
}

inline const scalar& Foam::diagStencil::center() const
{
    return this->v_;
}

inline scalar& Foam::diagStencil::center()
{
    return this->v_;
}

inline const scalar& diagStencil::operator[]
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= 1)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_;
}

inline scalar& diagStencil::operator[]
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= 1)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_;
}

inline void diagStencil::operator+=
(
    const diagStencil& ss
)
{
    this->v_ += ss.v_;
}

inline void diagStencil::operator-=
(
    const diagStencil& ss
)
{
    this->v_ -= ss.v_;
}

inline void diagStencil::operator=(const scalar s)
{
    this->v_ = s;
}

inline void diagStencil::operator=(const Foam::zero)
{
    this->v_ = Zero;
}

inline void diagStencil::operator*=
(
    const scalar s
)
{
    this->v_ *= s;
}

inline void diagStencil::operator/=
(
    const scalar s
)
{
    this->v_ /= s;
}

inline diagStencil operator+
(
    const diagStencil& st
)
{
    return st;
}

inline diagStencil operator-
(
    const diagStencil& st
)
{
    return diagStencil
    (
        -st[diagStencil::C]
    );
}

inline diagStencil operator+
(
    const diagStencil& ss1,
    const diagStencil& ss2
)
{
    return diagStencil
    (
        ss1[diagStencil::C]+ss2[diagStencil::C]
    );
}

inline diagStencil operator-
(
    const diagStencil& ss1,
    const diagStencil& ss2
)
{
    return diagStencil
    (
        ss1[diagStencil::C]-ss2[diagStencil::C]
    );
}

inline diagStencil operator*
(
    scalar s,
    const diagStencil& ss
)
{
    return diagStencil
    (
        s*ss[diagStencil::C]
    );
}

inline diagStencil operator*
(
    const diagStencil& ss,
    scalar s
)
{
    return s*ss;
}

inline diagStencil operator/
(
    const diagStencil& ss,
    scalar s
)
{
    return diagStencil
    (
        ss[diagStencil::C]/s
    );
}

inline bool operator==
(
    const diagStencil& st1,
    const diagStencil& st2
)
{
    return st1[diagStencil::C] == st2[diagStencil::C];
}

inline bool operator!=
(
    const diagStencil& st1,
    const diagStencil& st2
)
{
    return st1[diagStencil::C] != st2[diagStencil::C];
}

inline scalar stencilSum(const diagStencil& st)
{
    return st[diagStencil::C];
}

inline scalar neighborSum(const diagStencil& st)
{
    return 0.0;
}

inline scalar neighborSkewSum(const diagStencil& st)
{
    return 0.0;
}

template<>
class scalarProduct<scalar, diagStencil>
{
public:

    typedef diagStencil type;
};

template<>
class scalarProduct<diagStencil, scalar>
{
public:

    typedef diagStencil type;
};

}
