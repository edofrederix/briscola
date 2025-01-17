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

}
