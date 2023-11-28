#include "partPoints.H"

namespace Foam
{

namespace briscola
{

partPoints::partPoints()
:
    vectorBlock()
{}

partPoints::~partPoints()
{}

void partPoints::clear()
{
    vectorBlock::clear();
}

void partPoints::setSizeFromCells(const labelVector& size)
{
    // On each face of the block a layer of ghost points, plus one additional
    // layer because these are points, not cells.

    vectorBlock::setSize(size + 3*unitXYZ);
}

}

}
