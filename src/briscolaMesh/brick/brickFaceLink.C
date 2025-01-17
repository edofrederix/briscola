#include "brickFaceLink.H"
#include "brickTopology.H"
#include "brickLink.H"

namespace Foam
{

namespace briscola
{

brickFaceLink::brickFaceLink
(
    const face& f0,
    const face& f1,
    const bool periodic
)
:
    f0_(f0),
    f1_(f1),
    periodic_(periodic)
{
    // Compute transformation tensor

    const labelVector offset0 = faceOffsets[f0.num()];
    const labelVector offset1 = faceOffsets[f1.num()];

    T_ = offsetRotation(offset1, -offset0);

    brick b2(f1.parentBrick());

    b2.transform(T_);

    // Rotate around face normal until edges align

    const label k(faceNumber(-offset0));

    while (f0.e0() != b2.f(k).e0())
    {
        T_ = (b2.rotate(1, f0.num()/2) & T_);
    }

    // Check number of cells on either side

    labelVector N0 = f0.parentBrick().N();
    labelVector N1 = f1.parentBrick().N();

    N0[f0.num()/2] = 1;
    N1[f1.num()/2] = 1;

    if (N0 != cmptMag(T_ & N1))
    {
        FatalErrorInFunction
            << "Inconsistent number of cells across bricks "
            << f0.parentBrick().num() << " and "
            << f1.parentBrick().num() << " at "
            << f0 << endl
            << abort(FatalError);
    }
}

brickFaceLink::brickFaceLink
(
    const face& f0,
    const face& f1,
    const labelList& P,
    const brickTopology& topo,
    const bool periodic
)
:
    f0_(f0),
    f1_(f1),
    periodic_(periodic)
{
    // Compute transformation tensor by tracing back the non-periodic brick path

    T_ = eye;

    for (label i = P.size()-2; i >= 0; i--)
    {
        const brickLinks& links = topo.links()[P[i]];
        const brickFaceLink& link = links.getFaceLink(P[i+1]);

        T_ = link.T() & T_;
    }
}

brickFaceLink::~brickFaceLink()
{}

brickLink brickFaceLink::link() const
{
    return brickLink(*this);
}

}

}
