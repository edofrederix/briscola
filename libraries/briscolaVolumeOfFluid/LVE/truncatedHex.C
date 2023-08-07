#include "truncatedHex.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

void truncatedHex::setTets()
{
    tets_.clear();
    tets_.setSize(5);

    for (int i = 0; i < 5; i++)
    {
        vectorList v(4);

        for (int j = 0; j < 4; j++)
            v[j] = v_[tetDecomp[i][j]];

        tets_.set(i, new truncatedTet(v,n_,C_));
    }
}

truncatedHex::truncatedHex
(
    const vertexVector& v,
    const vector& n,
    const scalar& C
)
:
    v_(v),
    n_(n),
    C_(C)
{
    setTets();
}

truncatedHex::truncatedHex(const truncatedHex& hx)
:
    v_(hx.v_),
    n_(hx.n_),
    C_(hx.C_),
    tets_(hx.tets_)
{}

truncatedHex::~truncatedHex()
{}

scalar truncatedHex::volume() const
{
    scalar V = 0;

    forAll(tets_, i)
        V += tets_[i].volume();

    return V;
}

}

}

}
