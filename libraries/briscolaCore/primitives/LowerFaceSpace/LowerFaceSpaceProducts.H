namespace Foam
{

// Sum

template<class Cmpt1, class Cmpt2>
class typeOfSum<Cmpt1,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<LowerFaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class typeOfSum<LowerFaceSpace<Cmpt1>,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename typeOfSum<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

// Outer products

template<class Cmpt1, class Cmpt2>
class outerProduct<Cmpt1,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<LowerFaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class outerProduct<LowerFaceSpace<Cmpt1>,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename outerProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

// Inner product

template<class Cmpt1, class Cmpt2>
class innerProduct<Cmpt1,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<LowerFaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class innerProduct<LowerFaceSpace<Cmpt1>,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename innerProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

// Cross product

template<class Cmpt1, class Cmpt2>
class crossProduct<Cmpt1,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<LowerFaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class crossProduct<LowerFaceSpace<Cmpt1>,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename crossProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

// Scalar product

template<class Cmpt1, class Cmpt2>
class scalarProduct<Cmpt1,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<LowerFaceSpace<Cmpt1>,Cmpt2>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

template<class Cmpt1, class Cmpt2>
class scalarProduct<LowerFaceSpace<Cmpt1>,LowerFaceSpace<Cmpt2>>
{
public:

    typedef typename scalarProduct<Cmpt1,Cmpt2>::type productType;
    typedef LowerFaceSpace<productType> type;
};

}
