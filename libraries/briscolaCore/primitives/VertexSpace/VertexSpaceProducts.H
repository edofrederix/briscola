namespace Foam
{

// Sum

template<class Cmpt1, class Cmpt2>
class typeOfSum<Cmpt1,VertexSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<VertexSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<VertexSpace<Cmpt1>,VertexSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

// Outer products

template<class Cmpt1, class Cmpt2>
class outerProduct<Cmpt1,VertexSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<VertexSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<VertexSpace<Cmpt1>,VertexSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

// Inner product

template<class Cmpt1, class Cmpt2>
class innerProduct<Cmpt1,VertexSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<VertexSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<VertexSpace<Cmpt1>,VertexSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

// Cross product

template<class Cmpt1, class Cmpt2>
class crossProduct<Cmpt1,VertexSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<VertexSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<VertexSpace<Cmpt1>,VertexSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

// Scalar product

template<class Cmpt1, class Cmpt2>
class scalarProduct<Cmpt1,VertexSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<VertexSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<VertexSpace<Cmpt1>,VertexSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef VertexSpace<productType> type;
};

}