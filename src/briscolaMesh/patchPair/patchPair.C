#include "patchPair.H"
#include "geometry.H"

namespace Foam
{

namespace briscola
{

void patchPair::checkConsistency() const
{
    if (patch0_.facePtrs().size() != patch1_.facePtrs().size())
    {
        FatalErrorInFunction
            << "Paired patches " << patch0_.name() << " and " << patch1_.name()
            << " have different numbers of faces in them." << endl
            << exit(FatalError);
    }

    forAll(patch0_.facePtrs(), i)
    {
        const face& face1 = *patch0_.facePtrs()[i];
        const face& face2 = *patch1_.facePtrs()[i];

        if (face1.N() != face2.N())
        {
            FatalErrorInFunction
                << "Face " << i << " of paired patches " << patch0_.name()
                << " and " << patch1_.name() << " do not have the same number"
                << " of grid cells." << endl
                << exit(FatalError);
        }
    }
}

patchPair::patchPair
(
    const geometry& g,
    const label num,
    const patch& patch0,
    const patch& patch1
)
:
    meshObject<geometry>(g, num),
    patch0_(patch0),
    patch1_(patch1)
{
    checkConsistency();
}

patchPair::patchPair(const patchPair& p)
:
    meshObject<geometry>(p),
    patch0_(p.patch0()),
    patch1_(p.patch1())
{
    checkConsistency();
}

patchPair::~patchPair()
{}

}

}
