#include "decomposition.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(decomposition, 0);
defineRunTimeSelectionTable(decomposition, dictionary);

void decomposition::updateGlobalData()
{
    brickNumPerProc_.setSize(Pstream::nProcs());
    brickPartPerProc_.setSize(Pstream::nProcs());

    brickNumPerProc_[Pstream::myProcNo()] = myBrickNum_;
    brickPartPerProc_[Pstream::myProcNo()] = myBrickPart_;

    Pstream::gatherList(brickNumPerProc_);
    Pstream::gatherList(brickPartPerProc_);

    Pstream::scatterList(brickNumPerProc_);
    Pstream::scatterList(brickPartPerProc_);

    procMapPerBrick_.setSize(msh_.bricks().size());
    partSizePerBrick_.setSize(msh_.bricks().size());

    forAll(msh_.bricks(), bricki)
    {
        procMapPerBrick_.set(bricki, new labelBlock(decompPerBrick_[bricki]));
        procMapPerBrick_[bricki] = -1;

        partSizePerBrick_[bricki] =
            cmptDivide(msh_.bricks()[bricki].N(), decompPerBrick_[bricki]);
    }

    forAll(brickPartPerProc_, proci)
    {
        const label brickNum = brickNumPerProc_[proci];
        const labelVector brickPart = brickPartPerProc_[proci];

        procMapPerBrick_[brickNum](brickPart) = proci;
    }

    forAll(procMapPerBrick_, bricki)
    forAllBlock(procMapPerBrick_[bricki], i, j, k)
    {
        if (procMapPerBrick_[bricki](i,j,k) == -1)
        {
            FatalErrorInFunction
                << "Part index (" << i << "," << j << "," << k << ") is not "
                << "assigned to a processor" << endl
                << abort(FatalError);
        }
    }

}

decomposition::decomposition(mesh& msh)
:
    msh_(msh),
    dict_(msh.dict().subDict("decomposition")),
    decompPerBrick_(),
    myBrickNum_(),
    myBrickPart_(),
    brickNumPerProc_(),
    brickPartPerProc_(),
    procMapPerBrick_()
{}

decomposition::decomposition
(
    const decomposition& d
)
:
    msh_(d.msh_),
    dict_(d.dict_),
    decompPerBrick_(d.decompPerBrick_),
    myBrickNum_(d.myBrickNum_),
    myBrickPart_(d.myBrickPart_),
    brickNumPerProc_(d.brickNumPerProc_),
    brickPartPerProc_(d.brickPartPerProc_),
    procMapPerBrick_(d.procMapPerBrick_)
{}

decomposition::~decomposition()
{}

autoPtr<decomposition> decomposition::New(mesh& msh)
{
    const word decompType
    (
        msh.dict().subDict("decomposition").lookup("type")
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

    return autoPtr<decomposition>(cstrIter()(msh));
}

const brick& decomposition::myBrick() const
{
    return msh_.bricks()[myBrickNum_];
}


labelVector decomposition::myBrickN() const
{
    return myBrick().N();
}

}

}