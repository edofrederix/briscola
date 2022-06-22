#include "lineEdge.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace edges
{

defineTypeNameAndDebug(lineEdge, 0);
addToRunTimeSelectionTable(edge, lineEdge, dictionary);

lineEdge::lineEdge
(
    const face& f,
    const label num,
    const labelBlock& vertexNums,
    const labelVector& N
)
:
    edge(f, num, vertexNums, N)
{}

}

}

}
