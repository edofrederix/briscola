#include "tagAlgorithm.H"
#include "addToRunTimeSelectionTable.H"
#include "vof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(tagAlgorithm, 0);
defineRunTimeSelectionTable(tagAlgorithm, dictionary);

void tagAlgorithm::removeSmallParticles()
{
    colocatedLabelField& m = *this;

    labelList sizes(n_, Zero);

    const colocatedScalarField& cv =
        fvMsh_.template metrics<colocated>().cellVolumes();

    forAllCells(m,i,j,k)
    {
        if (m(i,j,k))
        {
            sizes[m(i,j,k) - 1]++;
        }
    }

    reduce(sizes, ListOp<sumOp<label>>());

    forAllCells(m,i,j,k)
    {
        if (m(i,j,k) && sizes[m(i,j,k) - 1] < minSize_)
        {
            removedVol_ += alpha_(i,j,k) * cv(i,j,k);

            m(i,j,k) = Zero;
            alpha_(i,j,k) = Zero;
        }
    }
}

tagAlgorithm::tagAlgorithm
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
        IOobject::NO_WRITE
    ),
    fvMsh_(fvMsh),
    dict_(dict),
    alpha_(alpha),
    n_(-1),
    minSize_(dict.lookupOrDefault<label>("minSize", 27)),
    removedVol_(Zero)
{
    // Initialize the label field to zero

    static_cast<colocatedLabelField&>(*this) = Zero;
}

tagAlgorithm::tagAlgorithm(const tagAlgorithm& s)
:
    colocatedLabelField(s),
    fvMsh_(s.fvMsh_),
    dict_(s.dict_),
    alpha_(s.alpha_),
    n_(s.n_),
    minSize_(s.minSize_),
    removedVol_(s.removedVol_)
{}

tagAlgorithm::~tagAlgorithm()
{}

autoPtr<tagAlgorithm> tagAlgorithm::New
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    colocatedScalarField& alpha
)
{
    const word tagAlgorithmType
    (
        dict.isNull()
      ? word("twoPass")
      : dict.lookup("type")
    );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(tagAlgorithmType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown tag algorithm type " << tagAlgorithmType
            << ". Valid tag algorithms are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<tagAlgorithm>
    (
        cstrIter()(fvMsh, dict, alpha)
    );
}

void tagAlgorithm::correct()
{
    binaryTag();

    this->tag();

    removeSmallParticles();
}

void tagAlgorithm::binaryTag()
{
    colocatedLabelField& m = *this;

    m = Zero;

    forAllCells(m,i,j,k)
        m(i,j,k) = alpha_(i,j,k) > 1e-3;
}

}

}

}
