namespace Foam
{

// With scalar

// Sum

template<>
class typeOfSum<scalar,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class typeOfSum<symmStencil,scalar>
{
public:

    typedef symmStencil type;
};

// Outer products

template<>
class outerProduct<scalar,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class outerProduct<symmStencil,scalar>
{
public:

    typedef symmStencil type;
};

// Inner product

template<>
class innerProduct<scalar,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class innerProduct<symmStencil,scalar>
{
public:

    typedef symmStencil type;
};

// Cross product

template<>
class crossProduct<scalar,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class crossProduct<symmStencil,scalar>
{
public:

    typedef symmStencil type;
};

// Scalar product

template<>
class scalarProduct<scalar,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class scalarProduct<symmStencil,scalar>
{
public:

    typedef symmStencil type;
};

// With self

// Sum

template<>
class typeOfSum<symmStencil,symmStencil>
{
public:

    typedef symmStencil type;
};

// Outer products

template<>
class outerProduct<symmStencil,symmStencil>
{
public:

    typedef symmStencil type;
};

// Inner product

template<>
class innerProduct<symmStencil,symmStencil>
{
public:

    typedef symmStencil type;
};

// Cross product

template<>
class crossProduct<symmStencil,symmStencil>
{
public:

    typedef symmStencil type;
};

// Scalar product

template<>
class scalarProduct<symmStencil,symmStencil>
{
public:

    typedef symmStencil type;
};

// With stencil

// Sum

template<>
class typeOfSum<stencil,symmStencil>
{
public:

    typedef stencil type;
};

template<>
class typeOfSum<symmStencil,stencil>
{
public:

    typedef stencil type;
};

// Outer products

template<>
class outerProduct<stencil,symmStencil>
{
public:

    typedef stencil type;
};

template<>
class outerProduct<symmStencil,stencil>
{
public:

    typedef stencil type;
};

// Inner product

template<>
class innerProduct<stencil,symmStencil>
{
public:

    typedef stencil type;
};

template<>
class innerProduct<symmStencil,stencil>
{
public:

    typedef stencil type;
};

// Cross product

template<>
class crossProduct<stencil,symmStencil>
{
public:

    typedef stencil type;
};

template<>
class crossProduct<symmStencil,stencil>
{
public:

    typedef symmStencil type;
};

// Scalar product

template<>
class scalarProduct<stencil,symmStencil>
{
public:

    typedef stencil type;
};

template<>
class scalarProduct<symmStencil,stencil>
{
public:

    typedef stencil type;
};

}
