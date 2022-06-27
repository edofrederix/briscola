namespace Foam
{

// Sum

template<class Cmpt1, class Cmpt2>
class typeOfSum<Cmpt1,FaceSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<FaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<FaceSpace<Cmpt1>,FaceSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

// Outer products

template<class Cmpt1, class Cmpt2>
class outerProduct<Cmpt1,FaceSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<FaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<FaceSpace<Cmpt1>,FaceSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

// Inner product

template<class Cmpt1, class Cmpt2>
class innerProduct<Cmpt1,FaceSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<FaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<FaceSpace<Cmpt1>,FaceSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

// Cross product

template<class Cmpt1, class Cmpt2>
class crossProduct<Cmpt1,FaceSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<FaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<FaceSpace<Cmpt1>,FaceSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

// Scalar product

template<class Cmpt1, class Cmpt2>
class scalarProduct<Cmpt1,FaceSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<FaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<FaceSpace<Cmpt1>,FaceSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef FaceSpace<productType> type;
};

}