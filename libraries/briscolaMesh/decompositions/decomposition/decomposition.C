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

    brickNumPerProc_[Pstream::myProcNo()] = myBrickNum();
    brickPartPerProc_[Pstream::myProcNo()] = myBrickPart();

    Pstream::gatherList(brickNumPerProc_);
    Pstream::gatherList(brickPartPerProc_);

    Pstream::scatterList(brickNumPerProc_);
    Pstream::scatterList(brickPartPerProc_);

    procMapPerBrick_.setSize(msh_.bricks().size());
    partSizePerBrick_.setSize(msh_.bricks().size());

    forAll(msh_.bricks(), bricki)
    {
        procMapPerBrick_.set
        (
            bricki,
            new labelBlock(decompPerBrick()[bricki], -1)
        );

        partSizePerBrick_[bricki] =
            cmptDivide(msh_.bricks()[bricki].N(), decompPerBrick()[bricki]);
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
        const labelBlock& brickMap = msh_.topology().map();

        // Get the number of processors per brick in each direction. For
        // structured brick topologies, the bricks are aligned with the local
        // coordinate system of the first brick. Also the brick map is aligned
        // with this coordinate system.

        labelList Nx(brickMap.l(), 0);
        labelList Ny(brickMap.m(), 0);
        labelList Nz(brickMap.n(), 0);

        forAllBlock(brickMap, i, j, k)
        {
            if (Nx[i] == 0 && brickMap(i,j,k) > -1)
            {
                Nx[i] = procMapPerBrick_[brickMap(i,j,k)].l();
            }

            if (Ny[j] == 0 && brickMap(i,j,k) > -1)
            {
                Ny[j] = procMapPerBrick_[brickMap(i,j,k)].m();
            }

            if (Nz[k] == 0 && brickMap(i,j,k) > -1)
            {
                Nz[k] = procMapPerBrick_[brickMap(i,j,k)].n();
            }
        }

        // Compute processor map, in a local coordinate system aligned with the
        // bricks and brick map.

        labelBlock map(sum(Nx), sum(Ny), sum(Nz), -1);

        labelVector cursor(zeroXYZ);

        for (int i = 0; i < brickMap.l(); i++)
        {
            cursor.y() = 0;

            for (int j = 0; j < brickMap.m(); j++)
            {
                cursor.z() = 0;

                for (int k = 0; k < brickMap.n(); k++)
                {
                    for (int ii = 0; ii < Nx[i]; ii++)
                    for (int jj = 0; jj < Ny[j]; jj++)
                    for (int kk = 0; kk < Nz[k]; kk++)
                    {
                        labelVector ijk(ii,jj,kk);

                        map(cursor+ijk) =
                            procMapPerBrick_[brickMap(i,j,k)](ijk);
                    }

                    cursor.z() += Nz[k];
                }

                cursor.y() += Ny[j];
            }

            cursor.x() += Nx[i];
        }

        // Store

        map_.setData(map);
    }
}

decomposition::decomposition(mesh& msh)
:
    msh_(msh),
    dict_(msh.dict().subDict("decomposition")),
    brickNumPerProc_(),
    brickPartPerProc_(),
    partSizePerProc_(),
    procMapPerBrick_(),
    map_()
{}

decomposition::decomposition
(
    const decomposition& d
)
:
    msh_(d.msh_),
    dict_(d.dict_),
    brickNumPerProc_(d.brickNumPerProc_),
    brickPartPerProc_(d.brickPartPerProc_),
    partSizePerProc_(d.partSizePerProc_),
    procMapPerBrick_(d.procMapPerBrick_),
    partSizePerBrick_(d.partSizePerBrick_),
    map_(d.map_)
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
    return msh_.bricks()[myBrickNum()];
}


labelVector decomposition::myBrickN() const
{
    return myBrick().N();
}

}

}