#include "decomposition.H"
#include "level.H"
#include "mesh.H"
#include "PstreamGlobals.H"

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

decomposition::decomposition(const decomposition& d)
:
    lvl_(d.lvl_),
    dict_(d.dict_),
    procBrickNums_(d.procBrickNums_),
    procBrickDecomps_(d.procBrickDecomps_),
    procBrickParts_(d.procBrickParts_),
    brickDecomps_(d.brickDecomps_),
    brickProcMaps_(d.brickProcMaps_),
    aggProcMap_(d.aggProcMap_),
    myAggPart_(d.myAggPart_),
    myAggProcNo_(d.myAggProcNo_),
    members_(d.members_),
    mapPtr_()
{
    if (d.mapPtr_.valid())
        mapPtr_.reset(new decompositionMap(d.map()));
}

decomposition::decomposition(const decomposition& d, const level& lvl)
:
    lvl_(lvl),
    dict_(d.dict_),
    procBrickNums_(d.procBrickNums_),
    procBrickDecomps_(d.procBrickDecomps_),
    procBrickParts_(d.procBrickParts_),
    brickDecomps_(d.brickDecomps_),
    brickProcMaps_(d.brickProcMaps_),
    aggProcMap_(d.aggProcMap_),
    myAggPart_(d.myAggPart_),
    myAggProcNo_(d.myAggProcNo_),
    members_(d.members_),
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
        if (procBrickParts_[proc] != -unitXYZ)
            brickProcMaps_[procBrickNums_[proc]](procBrickParts_[proc]) = proc;

    forAll(brickProcMaps_, b)
        forAllBlock(brickProcMaps_[b], i, j, k)
            if (brickProcMaps_[b](i,j,k) < 0)
                FatalErrorInFunction
                    << "Part " << labelVector(i,j,k) << " is not assigned "
                    << "to a processor" << endl << abort(FatalError);

    // Setup agglomerate data if needed, which is the case if the parent's brick
    // decomposition is not equal to the current brick decomposition

    if
    (
        lvl_.hasParent()
     && lvl_.parent().decomp().member()
     && lvl_.parent().decomp().myBrickDecomp() != myBrickDecomp()
    )
    {
        const decomposition& parent = lvl_.parent().decomp();

        const labelVector R =
            cmptDivide(parent.myBrickDecomp(), myBrickDecomp());

        labelVector start(parent.myBrickPart());

        for (int d = 0; d < 3; d++)
            start[d] = start[d] - (start[d] % R[d]);

        aggProcMap_.setSize(R);

        const labelBlock& parentMap =
            parent.brickProcMaps()[parent.myBrickNum()];

        // Slice the parent map over the agglomerate and store

        label proc = 0;
        for (label i = start.x(); i < start.x() + R.x(); i++)
        for (label j = start.y(); j < start.y() + R.y(); j++)
        for (label k = start.z(); k < start.z() + R.z(); k++)
        {
            labelVector ijk(i,j,k);

            aggProcMap_(ijk-start) = parentMap(ijk);

            if (parentMap(ijk) == Pstream::myProcNo())
            {
                myAggPart_ = ijk-start;
                myAggProcNo_ = proc;
            }

            proc++;
        }
    }

    // Set the decomposition members

    members_.resize(Pstream::nProcs());

    label proc = 0;
    forAll(brickProcMaps_, b)
        forAllBlock(brickProcMaps_[b], i, j, k)
            if (brickProcMaps_[b](i,j,k) > -1)
                members_[proc++] = brickProcMaps_[b](i,j,k);

    members_.resize(proc);

    // Set the global processor map if the brick topology is structured

    if (lvl_.msh().topology().structured())
        mapPtr_.reset(new decompositionMap(*this));
}

bool decomposition::coarsenable() const
{
    labelVector R = myBrickDecomp();

    if (R == unitXYZ)
    {
        return false;
    }
    else
    {
        for (int d = 0; d < 3; d++)
            R[d] /= Foam::max(R[d]/2, 1);

        const labelVector N(cmptMultiply(lvl_.N(), R));

        return coarsen(N) != N;
    }
}

bool decomposition::agglomerated() const
{
    return
        R() != unitXYZ
     && lvl_.hasParent()
     && lvl_.parent().decomp().member();
}

bool decomposition::aggParent() const
{
    return
        lvl_.hasChild()
     && lvl_.child().decomp().agglomerated();
}

}

}
