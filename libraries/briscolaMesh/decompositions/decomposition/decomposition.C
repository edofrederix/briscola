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

    // Set the part size per processor list

    partSizePerProc_.setSize(Pstream::nProcs());
    partSizePerProc_[Pstream::myProcNo()] = myPartN();

    Pstream::gatherList(partSizePerProc_);
    Pstream::scatterList(partSizePerProc_);

    // Set the global processor map if the brick topology is structured

    if (msh_.topology().structured())
    {
        const labelBlock& map = msh_.topology().map();

        // Get the number of processors per brick in each direction

        labelList Nx(map.l(), 0);
        labelList Ny(map.m(), 0);
        labelList Nz(map.n(), 0);

        forAllBlock(map, i, j, k)
        {
            if (Nx[i] == 0 && map(i,j,k) > -1)
            {
                Nx[i] = procMapPerBrick_[map(i,j,k)].l();
            }

            if (Ny[j] == 0 && map(i,j,k) > -1)
            {
                Ny[j] = procMapPerBrick_[map(i,j,k)].m();
            }

            if (Nz[k] == 0 && map(i,j,k) > -1)
            {
                Nz[k] = procMapPerBrick_[map(i,j,k)].n();
            }
        }

        // Set global processor map

        procMap_.setSize(sum(Nx), sum(Ny), sum(Nz));
        procMap_ = -1;

        labelVector cursor(zeroXYZ);

        for (int i = 0; i < map.l(); i++)
        {
            cursor.y() = 0;

            for (int j = 0; j < map.m(); j++)
            {
                cursor.z() = 0;

                for (int k = 0; k < map.n(); k++)
                {
                    for (int ii = 0; ii < Nx[i]; ii++)
                    for (int jj = 0; jj < Ny[j]; jj++)
                    for (int kk = 0; kk < Nz[k]; kk++)
                    {
                        labelVector ijk(ii,jj,kk);

                        procMap_(cursor+ijk) =
                            procMapPerBrick_[map(i,j,k)](ijk);
                    }

                    cursor.z() += Nz[k];
                }

                cursor.y() += Ny[j];
            }

            cursor.x() += Nx[i];
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
    partSizePerProc_(),
    procMapPerBrick_(),
    procMap_()
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
    partSizePerProc_(d.partSizePerProc_),
    procMapPerBrick_(d.procMapPerBrick_),
    partSizePerBrick_(d.partSizePerBrick_),
    procMap_(d.procMap_)
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