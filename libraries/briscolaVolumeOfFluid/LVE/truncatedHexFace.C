#include "truncatedHexFace.H"
#include "geometryFun.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

truncatedHexFace::truncatedHexFace(const vectorList& v)
:
    v_(v)
{}

truncatedHexFace::truncatedHexFace(const truncatedHexFace& hxf)
:
    v_(hxf.v_)
{}

void truncatedHexFace::operator=(const truncatedHexFace& hxf)
{
    v_ = hxf.v_;
}

truncatedHexFace::~truncatedHexFace()
{}

vector truncatedHexFace::vectorArea() const
{
    if (v_.size() == 3)
    {
        return triAreaVector(v_[0], v_[1], v_[2]);
    }
    else if (v_.size() == 4)
    {
        return
        (
            triAreaVector(v_[0], v_[1], v_[2])
          + triAreaVector(v_[2], v_[3], v_[0])
          + triAreaVector(v_[0], v_[1], v_[3])
          + triAreaVector(v_[1], v_[2], v_[3])
        ) / 2.0;
    }
    else if (v_.size() == 5)
    {
        return
        (
            // Cross triangles
            triAreaVector(v_[0], v_[2], v_[4])
          + triAreaVector(v_[0], v_[1], v_[3])
          + triAreaVector(v_[1], v_[2], v_[4])
          + triAreaVector(v_[2], v_[3], v_[0])
          + triAreaVector(v_[3], v_[4], v_[1])

            // Sides appear twice
          + 2.0*triAreaVector(v_[0], v_[1], v_[2])
          + 2.0*triAreaVector(v_[2], v_[3], v_[4])
          + 2.0*triAreaVector(v_[1], v_[2], v_[3])
          + 2.0*triAreaVector(v_[3], v_[4], v_[0])
          + 2.0*triAreaVector(v_[0], v_[1], v_[4])
        ) / 5.0;
    }
    else
    {
        FatalErrorInFunction
            << "Can only compute vector area for 3, 4 or 5 vertices"
            << endl << abort(FatalError);

        return vector::zero;
    }
}

vector truncatedHexFace::normal() const
{
    vector a(vectorArea());

    return a/Foam::mag(a);
}

scalar truncatedHexFace::area() const
{
    return Foam::mag(vectorArea());
}

}

}

}
