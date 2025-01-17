namespace Foam
{

// With scalar

// Sum

template<>
class typeOfSum<scalar,diagStencil>
{
public:

    typedef diagStencil type;
};

template<>
class typeOfSum<diagStencil,scalar>
{
public:

    typedef diagStencil type;
};

// Outer products

template<>
class outerProduct<scalar,diagStencil>
{
public:

    typedef diagStencil type;
};

template<>
class outerProduct<diagStencil,scalar>
{
public:

    typedef diagStencil type;
};

// Inner product

template<>
class innerProduct<scalar,diagStencil>
{
public:

    typedef diagStencil type;
};

template<>
class innerProduct<diagStencil,scalar>
{
public:

    typedef diagStencil type;
};

// Cross product

template<>
class crossProduct<scalar,diagStencil>
{
public:

    typedef diagStencil type;
};

template<>
class crossProduct<diagStencil,scalar>
{
public:

    typedef diagStencil type;
};

// Scalar product

template<>
class scalarProduct<scalar,diagStencil>
{
public:

    typedef diagStencil type;
};

template<>
class scalarProduct<diagStencil,scalar>
{
public:

    typedef diagStencil type;
};

// With self

// Sum

template<>
class typeOfSum<diagStencil,diagStencil>
{
public:

    typedef diagStencil type;
};

// Outer products

template<>
class outerProduct<diagStencil,diagStencil>
{
public:

    typedef diagStencil type;
};

// Inner product

template<>
class innerProduct<diagStencil,diagStencil>
{
public:

    typedef diagStencil type;
};

// Cross product

template<>
class crossProduct<diagStencil,diagStencil>
{
public:

    typedef diagStencil type;
};

// Scalar product

template<>
class scalarProduct<diagStencil,diagStencil>
{
public:

    typedef diagStencil type;
};

// With stencil

// Sum

template<>
class typeOfSum<stencil,diagStencil>
{
public:

    typedef stencil type;
};

template<>
class typeOfSum<diagStencil,stencil>
{
public:

    typedef stencil type;
};

// Outer products

template<>
class outerProduct<stencil,diagStencil>
{
public:

    typedef stencil type;
};

template<>
class outerProduct<diagStencil,stencil>
{
public:

    typedef stencil type;
};

// Inner product

template<>
class innerProduct<stencil,diagStencil>
{
public:

    typedef stencil type;
};

template<>
class innerProduct<diagStencil,stencil>
{
public:

    typedef stencil type;
};

// Cross product

template<>
class crossProduct<stencil,diagStencil>
{
public:

    typedef stencil type;
};

template<>
class crossProduct<diagStencil,stencil>
{
public:

    typedef diagStencil type;
};

// Scalar product

template<>
class scalarProduct<stencil,diagStencil>
{
public:

    typedef stencil type;
};

template<>
class scalarProduct<diagStencil,stencil>
{
public:

    typedef stencil type;
};

}
