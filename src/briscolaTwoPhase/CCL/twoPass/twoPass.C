#include "twoPass.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(twoPass, 0);
addToRunTimeSelectionTable(CCL, twoPass, dictionary);

twoPass::twoPass
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    CCL(fvMsh, dict, alpha)
{}

twoPass::twoPass(const twoPass& s)
:
    CCL(s)
{}

twoPass::~twoPass()
{}

void twoPass::tag()
{
    colocatedLabelField& m = *this;
}

}

}

}
