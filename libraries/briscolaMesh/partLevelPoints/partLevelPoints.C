#include "partLevelPoints.H"

namespace Foam
{

namespace briscola
{

partLevelPoints::partLevelPoints()
:
    vectorBlock()
{}

partLevelPoints::~partLevelPoints()
{}

void partLevelPoints::clear()
{
    vectorBlock::clear();
}

void partLevelPoints::setSizeFromCells(const labelVector& size)
{
    // On each face of the block a layer of ghost points, plus one additional
    // layer because these are points, not cells.

    vectorBlock::setSize(size + 3*unitXYZ);
}

}

}
