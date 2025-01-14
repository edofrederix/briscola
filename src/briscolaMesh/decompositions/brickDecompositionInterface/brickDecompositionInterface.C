#include "brickDecompositionInterface.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

void brickDecompositionInterface::setSlices()
{
    const labelVector offset = link_.offset();
    const labelTensor T = link_.T();

    // Processor topology of both bricks

    labelBlock procNums0(msh_.decomp().procMapPerBrick()[b0_.num()]);
    labelBlock procNums1(msh_.decomp().procMapPerBrick()[b1_.num()]);

    // Part sizes of both bricks

    labelVector partN0(msh_.decomp().partSizePerBrick()[b0_.num()]);
    labelVector partN1(msh_.decomp().partSizePerBrick()[b1_.num()]);

    // Sizes of both bricks

    labelVector Nb0 = b0_.N();
    labelVector Nb1 = b1_.N();

    // Transform data of the second brick to oppose the first brick

    procNums1.transform(T);
    partN1 = cmptMag(T & partN1);
    Nb1 = cmptMag(T & Nb1);

    // Slice the processor number blocks at the offset

    for (label dir = 0; dir < 3; dir++)
    if (offset[dir] != 0)
    {
        procNums0 = procNums0.slice(-(1+offset[dir])/2, dir);
        procNums1 = procNums1.slice(-(1-offset[dir])/2, dir);
    }

    if (procNums0.shape() != procNums1.shape())
    {
        FatalError
            << "Mismatch between the decompositions of brick "
            << link_.b0().num() << " and " << link_.b1().num()
            << " at offset " << offset << endl;
        FatalError.exit();
    }

    map_.setSize(procNums0.shape());
    slices_.setSize(map_.size());

    // Add to slices

    label l = 0;

    forAllBlock(map_, i, j, k)
    {
        const labelVector ijk(i,j,k);

        // Start and end position of this slice, relative to the brick

        const labelVector start =
            cmptMultiply(ijk, partN0)
          + cmptMultiply(cmptMax(offset,zeroXYZ), Nb0-unitXYZ);

        const labelVector end
        (
            start.x() + (offset.x() == 0 ? partN0.x() : 1),
            start.y() + (offset.y() == 0 ? partN0.y() : 1),
            start.z() + (offset.z() == 0 ? partN0.z() : 1)
        );

        const labelVector N(end-start);

        labelVector start0 = start;
        labelVector start1 = start;

        labelVector N0 = N;
        labelVector N1 = N;

        // Rotate and shift the slice in the second brick back to the second
        // brick's original orientation

        const labelTensor S(shift(T.T()));

        start1 = (T.T() & start1) + (S & (b1_.N() - unitXYZ));
        N1 = (T.T() & N1);

        // Assure that slices are directed in the positive direction

        start0 = cmptMin(start0, start0 + N0 + unitXYZ);
        start1 = cmptMin(start1, start1 + N1 + unitXYZ);

        N0 = cmptMag(N0);
        N1 = cmptMag(N1);

        const label procNum0(procNums0(ijk));
        const label procNum1(procNums1(ijk));

        // Store

        map_(ijk) = l;

        slices_.set
        (
            l,
            new brickDecompositionSlice
            (
                procNum0,
                procNum1,
                start0,
                start1,
                N0,
                N1
            )
        );

        l++;
    }
}

brickDecompositionInterface::brickDecompositionInterface
(
    const brickLink& link,
    const mesh& msh
)
:
    link_(link),
    b0_(link.b0()),
    b1_(link.b1()),
    msh_(msh)
{
    if (link.offset() != zeroXYZ)
    {
        setSlices();
    }
}

brickDecompositionInterface::~brickDecompositionInterface()
{}

}

}
