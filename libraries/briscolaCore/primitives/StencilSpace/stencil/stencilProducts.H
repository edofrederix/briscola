namespace Foam
{

// With scalar

// Sum

template<>
class typeOfSum<scalar,stencil>
{
public:

    typedef stencil type;
};

template<>
class typeOfSum<stencil,scalar>
{
public:

    typedef stencil type;
};

// Outer products

template<>
class outerProduct<scalar,stencil>
{
public:

    typedef stencil type;
};

template<>
class outerProduct<stencil,scalar>
{
public:

    typedef stencil type;
};

// Inner product

template<>
class innerProduct<scalar,stencil>
{
public:

    typedef stencil type;
};

template<>
class innerProduct<stencil,scalar>
{
public:

    typedef stencil type;
};

// Cross product

template<>
class crossProduct<scalar,stencil>
{
public:

    typedef stencil type;
};

template<>
class crossProduct<stencil,scalar>
{
public:

    typedef stencil type;
};

// Scalar product

template<>
class scalarProduct<scalar,stencil>
{
public:

    typedef stencil type;
};

template<>
class scalarProduct<stencil,scalar>
{
public:

    typedef stencil type;
};

// With self

// Sum

template<>
class typeOfSum<stencil,stencil>
{
public:

    typedef stencil type;
};

// Outer products

template<>
class outerProduct<stencil,stencil>
{
public:

    typedef stencil type;
};

// Inner product

template<>
class innerProduct<stencil,stencil>
{
public:

    typedef stencil type;
};

// Cross product

template<>
class crossProduct<stencil,stencil>
{
public:

    typedef stencil type;
};

// Scalar product

template<>
class scalarProduct<stencil,stencil>
{
public:

    typedef stencil type;
};

// With symmetric stencil

// Sum

template<>
class typeOfSum<stencil,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class typeOfSum<symmStencil,stencil>
{
public:

    typedef symmStencil type;
};

// Outer products

template<>
class outerProduct<stencil,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class outerProduct<symmStencil,stencil>
{
public:

    typedef symmStencil type;
};

// Inner product

template<>
class innerProduct<stencil,symmStencil>
{
public:

    typedef symmStencil type;
};

template<>
class innerProduct<symmStencil,stencil>
{
public:

    typedef symmStencil type;
};

// Cross product

template<>
class crossProduct<stencil,symmStencil>
{
public:

    typedef symmStencil type;
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

    typedef symmStencil type;
};

template<>
class scalarProduct<symmStencil,stencil>
{
public:

    typedef symmStencil type;
};

}
