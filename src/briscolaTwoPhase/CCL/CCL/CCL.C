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

void CCL::removeSmallComponents()
{
    colocatedLabelField& m = *this;

    labelList sizes(n_, Zero);

    forAllCells(m,i,j,k)
    {
        if (m(i,j,k))
        {
            sizes[m(i,j,k) - 1]++;
        }
    }

    forAll(sizes, i)
    {
        reduce(sizes[i], sumOp<label>());
    }

    forAllCells(m,i,j,k)
    {
        if (m(i,j,k) && sizes[m(i,j,k) - 1] < minSize_)
        {
            m(i,j,k) = Zero;
            alpha_(i,j,k) = Zero;
        }
    }
}

CCL::CCL
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    colocatedScalarField& alpha
)
:
    colocatedLabelField
    (
        "m",
        fvMsh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        false,
        false
    ),
    fvMsh_(fvMsh),
    dict_(dict),
    alpha_(alpha),
    n_(-1),
    minSize_(dict.lookupOrDefault<label>("minSize", 27))
{
    // Initialize the label field to zero

    static_cast<colocatedLabelField&>(*this) = Zero;
}

CCL::CCL(const CCL& s)
:
    colocatedLabelField(s),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    alpha_(s.alpha_),
    n_(s.n_),
    minSize_(s.minSize_)
{}

CCL::~CCL()
{}

autoPtr<CCL> CCL::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    colocatedScalarField& alpha
)
{
    const word CCLType
    (
        dict.isNull()
      ? word("twoPass")
      : dict.lookup("type")
    );

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

void CCL::correct()
{
    binaryTag();

    this->tag();

    removeSmallComponents();
}

void CCL::binaryTag()
{
    colocatedLabelField& m = *this;

    m = Zero;

    forAllCells(m,l,d,i,j,k)
        m(l,d,i,j,k) = alpha_(l,d,i,j,k) > 1e-3;
}

}

}

}
