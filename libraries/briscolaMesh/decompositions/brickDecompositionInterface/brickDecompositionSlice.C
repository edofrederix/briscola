#include "brickDecompositionSlice.H"

namespace Foam
{

namespace briscola
{

brickDecompositionSlice::brickDecompositionSlice
(
    const label procNum0,
    const label procNum1,
    const labelVector& start0,
    const labelVector& start1,
    const labelVector& N0,
    const labelVector& N1
)
:
    procNum0_(procNum0),
    procNum1_(procNum1),
    start0_(start0),
    start1_(start1),
    N0_(N0),
    N1_(N1)
{}

brickDecompositionSlice::~brickDecompositionSlice()
{}

}

}
