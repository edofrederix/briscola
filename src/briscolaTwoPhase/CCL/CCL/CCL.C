#include "CCL.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(CCL, 0);
defineRunTimeSelectionTable(CCL, dictionary);

CCL::CCL
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
:
    colocatedLabelField
    (
        "m",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false
    ),
    fvMsh_(fvMsh),
    dict_(dict),
    alpha_(alpha)
{
    // Initialize the label field to zero

    static_cast<colocatedLabelField&>(*this) = Zero;
}

CCL::CCL(const CCL& s)
:
    colocatedLabelField(s),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    alpha_(s.alpha_)
{}

CCL::~CCL()
{}

autoPtr<CCL> CCL::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
)
{
    const word CCLType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(CCLType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown CCL algorithm type " << CCLType
            << ". Valid CCL algorithm are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<CCL>
    (
        cstrIter()(fvMsh, dict, alpha)
    );
}

void CCL::binaryTag()
{
    colocatedLabelField& m = *this;

    forAllCells(m,l,d,i,j,k)
        m(l,d,i,j,k) = alpha_(l,d,i,j,k) > vof::threshold;
}

}

}

}
