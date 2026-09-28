#include "lubricationBoundaryCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(lubricationBoundaryCondition, 0);

typedef boundaryCondition<vector,colocated> colocatedVectorBoundaryCondition;

addToRunTimeSelectionTable
(
    colocatedVectorBoundaryCondition,
    lubricationBoundaryCondition,
    dictionary
);

lubricationBoundaryCondition::lubricationBoundaryCondition
(
    const colocatedVectorLevel& level,
    const boundary& b
)
:
    NeumannBoundaryCondition<vector,colocated>(level, b, Zero),
    inverted_(this->dict().lookupOrDefault<Switch>("inverted", false))
{}

lubricationBoundaryCondition::lubricationBoundaryCondition
(
    const lubricationBoundaryCondition& bc
)
:
    NeumannBoundaryCondition<vector,colocated>(bc),
    inverted_(bc.inverted_)
{}

lubricationBoundaryCondition::lubricationBoundaryCondition
(
    const lubricationBoundaryCondition& bc,
    const colocatedVectorLevel& level
)
:
    NeumannBoundaryCondition<vector,colocated>(bc, level),
    inverted_(bc.inverted_)
{}

void lubricationBoundaryCondition::evaluate(const label)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label l(this->l_);
    const label fd = f/2;

    // By definition, the interface normal points into regions of high alpha.
    // The boundary condition is designed to keep alpha away from the boundary.
    // So, the normal should be set as the face normal, but pointing into the
    // domain. Face normals of upper boundaries point into the domain but those
    // of lower boundaries point outward and need to have their signs changed.

    const label sign = (f%1 ? -1 : +1)*(inverted_ ? -1 : +1);

    colocatedVectorDirection& normal = this->level_[0];

    const colocatedVectorFaceField& fn = this->faceNormals();

    const labelVector S(this->S(0));
    const labelVector E(this->E(0));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const labelVector upp(upperFaceNeighbor(ijk,f));

        // Set the internal value equal to the boundary face normal

        if (Foam::mag(normal(ijk)))
            normal(ijk) = sign*fn[fd](l,0,upp);
    }
}

}

}

}

