#include "decomposition.H"
#include "level.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(decomposition, 0);
defineRunTimeSelectionTable(decomposition, dictionary);

decomposition::decomposition(const level& lvl)
:
    lvl_(lvl),
    dict_(lvl.msh().dict().subDict("decomposition"))
{}

decomposition::decomposition
(
    const decomposition& d
)
:
    lvl_(d.lvl_),
    dict_(d.dict_),
    procBrickNums_(d.procBrickNums_),
    procBrickDecomps_(d.procBrickDecomps_),
    procBrickParts_(d.procBrickParts_),
    brickDecomps_(d.brickDecomps_),
    brickProcMaps_(d.brickProcMaps_),
    mapPtr_()
{
    if (d.mapPtr_.valid())
        mapPtr_.reset(new decompositionMap(d.map()));
}

decomposition::~decomposition()
{}

autoPtr<decomposition> decomposition::New(const level& lvl)
{
    const word decompType
    (
        lvl.msh().dict().subDict("decomposition").lookup("type")
    );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(decompType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown decomposition type " << decompType
            << ". Valid decomposition types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<decomposition> ptr(cstrIter()(lvl));

    ptr->init();

    return ptr;
}

void decomposition::init()
{
    // Collect processor values

    procBrickNums_.setSize(Pstream::nProcs());
    procBrickDecomps_.setSize(Pstream::nProcs());
    procBrickParts_.setSize(Pstream::nProcs());

    procBrickNums_[Pstream::myProcNo()] = myBrickNum();
    procBrickDecomps_[Pstream::myProcNo()] = myBrickDecomp();
    procBrickParts_[Pstream::myProcNo()] = myBrickPart();

    Pstream::gatherList(procBrickNums_);
    Pstream::gatherList(procBrickDecomps_);
    Pstream::gatherList(procBrickParts_);

    Pstream::scatterList(procBrickNums_);
    Pstream::scatterList(procBrickDecomps_);
    Pstream::scatterList(procBrickParts_);

    // Store brick decompositions

    brickDecomps_.setSize(lvl_.msh().bricks().size());

    forAll(procBrickNums_, proc)
        brickDecomps_[procBrickNums_[proc]] =
            procBrickDecomps_[proc];

    // Make brick processor map

    brickProcMaps_.setSize(lvl_.msh().bricks().size());

    forAll(brickDecomps_, b)
        brickProcMaps_.set(b, new labelBlock(brickDecomps_[b], -1));

    forAll(procBrickNums_, proc)
        brickProcMaps_[procBrickNums_[proc]](procBrickParts_[proc]) = proc;

    forAll(brickProcMaps_, b)
        forAllBlock(brickProcMaps_[b], i, j, k)
            if (brickProcMaps_[b](i,j,k) < 0)
                FatalErrorInFunction
                    << "Part " << labelVector(i,j,k) << " is not assigned "
                    << "to a processor" << endl << abort(FatalError);

    // Set the global processor map if the brick topology is structured

    if (lvl_.msh().topology().structured())
        mapPtr_.reset(new decompositionMap(*this));
}

}

}
