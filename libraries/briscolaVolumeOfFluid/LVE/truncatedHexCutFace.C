#include "truncatedHexCutFace.H"
#include "geometryFun.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

truncatedHexCutFace::truncatedHexCutFace(const vectorList& v)
:
    truncatedHexFace(v)
{}

truncatedHexCutFace::truncatedHexCutFace(const truncatedHexCutFace& hxf)
:
    truncatedHexFace(hxf)
{}

truncatedHexCutFace::~truncatedHexCutFace()
{}

vector truncatedHexCutFace::vectorArea() const
{
    if (v_.size() == 3)
    {
        return triAreaVector(v_[0], v_[1], v_[2]);
    }
    else if (v_.size() == 4)
    {
        return
            triAreaVector(v_[0], v_[1], v_[2])
          + triAreaVector(v_[2], v_[3], v_[0]);
    }
    else if (v_.size() == 5)
    {
        return
            triAreaVector(v_[0], v_[1], v_[2])
          + triAreaVector(v_[2], v_[3], v_[4])
          + triAreaVector(v_[0], v_[2], v_[4]);
    }
    else if (v_.size() == 6)
    {
        return
            triAreaVector(v_[0], v_[1], v_[2])
          + triAreaVector(v_[2], v_[3], v_[4])
          + triAreaVector(v_[4], v_[5], v_[0])
          + triAreaVector(v_[0], v_[2], v_[4]);
    }
    else
    {
        FatalErrorInFunction
            << "Can only compute vector area for 3, 4, 5 or 6 vertices"
            << endl << abort(FatalError);

        return vector::zero;
    }
}

}

}

}
