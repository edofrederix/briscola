#include "decompositionInterface.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

void decompositionInterface::setSlices()
{
    const labelVector offset = link_.offset();
    const labelTensor T = link_.T();

    // Processor topology of both bricks

    labelBlock procNums0 = lvl_.decomp().brickProcMaps()[b0_.num()];
    labelBlock procNums1 = lvl_.decomp().brickProcMaps()[b1_.num()];

    // Transform data of the second brick to oppose the first brick

    procNums1.transform(T);

    // Slice the processor number blocks at the offset

    for (label dir = 0; dir < 3; dir++)
    if (offset[dir] != 0)
    {
        procNums0 = procNums0.slice(-(1+offset[dir])/2, dir);
        procNums1 = procNums1.slice(-(1-offset[dir])/2, dir);
    }

    if (procNums0.shape() != procNums1.shape())
        FatalErrorInFunction
            << "Mismatch between the decompositions of brick "
            << link_.b0().num() << " and " << link_.b1().num()
            << " at offset " << offset << endl
            << abort(FatalError);

    map_.setSize(procNums0.shape());
    slices_.setSize(map_.size());

    // Add to slices

    label l = 0;

    forAllBlock(map_, i, j, k)
    {
        const labelVector ijk(i,j,k);

        map_(ijk) = l;

        slices_.set(l, new labelPair(procNums0(ijk), procNums1(ijk)));

        l++;
    }
}

decompositionInterface::decompositionInterface
(
    const brickLink& link,
    const level& lvl
)
:
    lvl_(lvl),
    link_(link),
    b0_(link.b0()),
    b1_(link.b1())
{
    if (link.offset() != zeroXYZ)
    {
        setSlices();
    }
}

}

}
