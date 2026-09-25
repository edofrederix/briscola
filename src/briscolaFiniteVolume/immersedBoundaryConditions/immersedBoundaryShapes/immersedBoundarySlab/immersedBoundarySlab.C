#include "immersedBoundarySlab.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

immersedBoundarySlab::immersedBoundarySlab
(
    const dictionary& dict,
    bool inverted
)
:
    immersedBoundaryShape(dict,inverted),
    center_(vector(dict.lookup("center"))),
    normal_(vector(dict.lookup("normal"))),
    thickness_(readScalar(dict.lookup("thickness")))
{
    if (mag(normal_) < SMALL)
        FatalErrorInFunction
            << "Slab normal must have a non-zero length."
            << exit(FatalError);

    normal_ /= mag(normal_);

    if (thickness_ <= 0.0)
        FatalErrorInFunction
            << "Slab thickness has to be a positive value."
            << exit(FatalError);
}

// Destructor

immersedBoundarySlab::~immersedBoundarySlab()
{}

bool immersedBoundarySlab::isInside(vector point) const
{
    // Distance to the closest face (negative if inside slab)
    const scalar dist = Foam::mag((point - center_) & normal_) - 0.5*thickness_;

    // Points exactly on the slab's surface are treated as ib
    return inverted_ ? dist > 0 : dist < 0;
}

scalar immersedBoundarySlab::wallDistance(vector c, vector nb) const
{
    // Return -1 if the center point is not a fluid point or if the neighboring
    // point is not inside the shape
    if (this->isInside(c) || !this->isInside(nb))
        return -1;

    const scalar sc  = ((c  - center_) & normal_);
    const scalar snb = ((nb - center_) & normal_);

    const scalar sf = Foam::sign(inverted_ ? snb : sc)*0.5*thickness_;
    const scalar f = Foam::max((sf - sc)/(snb - sc), scalar(0));

    return f*Foam::mag(nb - c);
}

scalar immersedBoundarySlab::wallNormalDistance(vector gc) const
{
    // Return -1 if the point is not inside the shape
    if (!this->isInside(gc))
        return -1;

    const scalar sg = ((gc - center_) & normal_);

    return inverted_
      ? Foam::mag(sg) - 0.5*thickness_
      : 0.5*thickness_ - Foam::mag(sg);
}

vector immersedBoundarySlab::mirrorPoint(vector gc) const
{
    if (this->wallNormalDistance(gc) <= 0)
        return gc;

    const scalar sg = ((gc - center_) & normal_);

    // Reflect in the face on the same side of the mid-plane as gc
    return gc + 2.0*(Foam::sign(sg)*0.5*thickness_ - sg)*normal_;
}

}

}

}
