#include "MittalDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>::
MittalDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>
    (
        mshField,
        ib,
        ib.ghostMask()
    ),
    exchangePoints_(this->IB_.mask().numberOfLevels()),
    exchanges_(this->IB_.mask().numberOfLevels()),
    boundaryValues_(this->dict().lookup("values"))
{
    // Check shape overlap
    if (this->IB_.shapeOverlap())
    {
        WarningInFunction
            << "Overlapping shapes identified."
            << " This may cause issues with Mittal IBM." << endl;
    }

    forAll(exchangePoints_, l)
    {
        exchangePoints_[l].setSize(MeshType::numberOfDirections);
        exchanges_[l].setSize(MeshType::numberOfDirections);
    }

    // Mesh
    const mesh& msh = this->fvMsh_.msh();

    // Set exchange points
    forAllCells(this->IB_.mirrorPoints(),l,d,i,j,k)
    {
        if (this->forcingPoints_(l,d,i,j,k))
        {
            vector mp = this->IB_.mirrorPoints()(l,d,i,j,k);

            if (this->IB_.isInside(mp))
            {
                WarningInFunction
                    << "Mirror point of ghost cell " << vector(i,j,k)
                    << " at d = "<< d << " located inside of "
                    << "immersed boundary."
                    << " This may cause issues with Mittal IBM." << endl;
            }

            // Colocated cell index of mp
            labelVector mpIndex = msh.findCell(mp, l);

            if (mpIndex == -unitXYZ)
            {
                exchangePoints_[l][d].append(mp);
            }
        }
    }

    // Set data point exchanges
    forAll(exchanges_, l)
    {
        forAll(exchanges_[l], d)
        {
            exchanges_[l].set
            (
                d,
                new pointDataExchange<MeshType>
                (
                    exchangePoints_[l][d],
                    this->fvMsh_,
                    l,
                    d
                )
            );
        }
    }
}

// Destructor

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>::
~MittalDirichletImmersedBoundaryCondition()
{}

template<class Type, class MeshType>
void MittalDirichletImmersedBoundaryCondition<Type,MeshType>
::correctJacobiPoints
(
    meshLevel<Type,MeshType>& x
) const
{
    scalar omega = this->omega_;

    // Mesh
    const fvMesh& fvMsh = x.fvMsh();
    const mesh& msh = fvMsh.msh();

    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh.metrics<MeshType>().cellCenters();

    label l = x.levelNum();

    forAll(x, d)
    {
        List<Type> exchangeData
        (
            move(exchanges_[l][d].dataFunc(x.mshField()))
        );

        scalar cursor = 0;

        forAllCells(x[d],i,j,k)
        {
            if (this->forcingPoints_(l,d,i,j,k))
            {
                // mirror point
                vector mp = this->IB_.mirrorPoints()(l,d,i,j,k);

                // Colocated cell index of mp
                labelVector mpIndex = msh.findCell(mp, l);

                Type mpValue = Zero;

                if (mpIndex == -unitXYZ)
                {
                    mpValue = exchangeData[cursor++];
                }
                else
                {
                    // Local coordinates of mp in colocated cell
                    vector mpLocalCoords = msh[l].points()
                        .cellCoordinates(mp, mpIndex, true);

                    if (mpLocalCoords == -vector::one)
                    {
                        FatalError
                            << "Interpolation error"
                            << " at direction " << d
                            << ". Mirror point: " << mp
                            << " and colocated cell index: "
                            << mpIndex
                            << endl;
                        FatalError.exit();
                    }

                    // Index of <MeshType> left-bottom-aft cell w.r.t. mp
                    labelVector mpLBA = mpIndex;

                    for (int dir = 0; dir < 3; dir++)
                    {
                        if
                        (
                               (!MeshType::padding[d][dir])
                            && (mpLocalCoords[dir] < 0.5)
                        )
                        {
                            mpLBA[dir] -= 1;
                        }
                    }

                    // Interpolation box
                    vertexVector interpPoints
                    (
                        CC[l][d](mpLBA),
                        CC[l][d](mpLBA+unitX),
                        CC[l][d](mpLBA+unitY),
                        CC[l][d](mpLBA+unitXY),
                        CC[l][d](mpLBA+unitZ),
                        CC[l][d](mpLBA+unitXZ),
                        CC[l][d](mpLBA+unitYZ),
                        CC[l][d](mpLBA+unitXYZ)
                    );

                    // Interpolation weights
                    const vector v
                    (
                        interpolationWeights(mp,interpPoints,true)
                    );

                    vertexScalar weights
                    (
                        (1-v.x())*(1-v.y())*(1-v.z()),
                        (  v.x())*(1-v.y())*(1-v.z()),
                        (1-v.x())*(  v.y())*(1-v.z()),
                        (  v.x())*(  v.y())*(1-v.z()),
                        (1-v.x())*(1-v.y())*(  v.z()),
                        (  v.x())*(1-v.y())*(  v.z()),
                        (1-v.x())*(  v.y())*(  v.z()),
                        (  v.x())*(  v.y())*(  v.z())
                    );

                    if (v != -vector::one)
                    {
                        mpValue =
                              weights.lba()*x[d](mpLBA)
                            + weights.rba()*x[d](mpLBA+unitX)
                            + weights.lta()*x[d](mpLBA+unitY)
                            + weights.rta()*x[d](mpLBA+unitXY)
                            + weights.lbf()*x[d](mpLBA+unitZ)
                            + weights.rbf()*x[d](mpLBA+unitXZ)
                            + weights.ltf()*x[d](mpLBA+unitYZ)
                            + weights.rtf()*x[d](mpLBA+unitXYZ);
                    }
                    else
                    {
                        FatalError
                            << "Interpolation error"
                            << " at direction " << d
                            << ". Mirror point: " << mp << nl
                            << "and interpolation points: "
                            << interpPoints << nl
                            << "Local colocated coordinates: "
                            << mpLocalCoords << nl
                            << "LBA cell index: " << mpLBA << nl
                            << "Colocated cell index" << mpIndex
                            << endl;
                        FatalError.exit();
                    }
                }

                // If the ghost cell is on the boundary or at the center of
                // a sphere or cylinder, set value to boundary value

                if (mp == CC(l,d,i,j,k))
                {
                    mpValue = boundaryValues_[d];
                }

                x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k)
                    + omega * (2.0*boundaryValues_[d] - mpValue);
            }
        }
    }
}

}

}

}
