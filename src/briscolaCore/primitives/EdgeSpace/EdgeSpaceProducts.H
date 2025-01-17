namespace Foam
{

// Sum

template<class Cmpt1, class Cmpt2>
class typeOfSum<Cmpt1,EdgeSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<EdgeSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<EdgeSpace<Cmpt1>,EdgeSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

// Outer products

template<class Cmpt1, class Cmpt2>
class outerProduct<Cmpt1,EdgeSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<EdgeSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<EdgeSpace<Cmpt1>,EdgeSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

// Inner product

template<class Cmpt1, class Cmpt2>
class innerProduct<Cmpt1,EdgeSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<EdgeSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<EdgeSpace<Cmpt1>,EdgeSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

// Cross product

template<class Cmpt1, class Cmpt2>
class crossProduct<Cmpt1,EdgeSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<EdgeSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<EdgeSpace<Cmpt1>,EdgeSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

// Scalar product

template<class Cmpt1, class Cmpt2>
class scalarProduct<Cmpt1,EdgeSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<EdgeSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<EdgeSpace<Cmpt1>,EdgeSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef EdgeSpace<productType> type;
};

}