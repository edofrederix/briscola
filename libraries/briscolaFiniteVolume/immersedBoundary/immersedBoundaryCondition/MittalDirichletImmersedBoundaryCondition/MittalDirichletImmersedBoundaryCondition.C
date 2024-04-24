#include "MittalDirichletImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>
::MittalDirichletImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& mshField,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(mshField,ib,true),
    exchangePoints_(this->IB_.mask().numberOfLevels()),
    boundaryValues_(this->dict().lookup("values"))
{
    forAll(exchangePoints_, l)
    {
        exchangePoints_[l].setSize(MeshType::numberOfDirections);
    }

    // Mesh
    const mesh& msh = this->fvMsh_.msh();

    forAllCells(this->IB_.mirrorPoints(),l,d,i,j,k)
    {
        if (this->IB_.ghostMask()(l,d,i,j,k))
        {
            vector mp = this->IB_.mirrorPoints()(l,d,i,j,k);

            // Colocated cell index of mp
            labelVector mpIndex = msh.findCell(mp, l);

            if (mpIndex == -unitXYZ)
            {
                exchangePoints_[l][d].append(mp);
            }
        }
    }
}

// Destructor

template<class Type, class MeshType>
MittalDirichletImmersedBoundaryCondition<Type,MeshType>
::~MittalDirichletImmersedBoundaryCondition()
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

    const scalar H = l == 0;

    forAll(x, d)
    {
        pointDataExchange<MeshType> exchange(exchangePoints_[l][d], fvMsh, l, d);

        List<Type> exchangeData(move(exchange(x.mshField())));

        scalar cursor = 0;

        forAllCells(x[d],i,j,k)
        {
            if (this->IB_.ghostMask()(l,d,i,j,k))
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
                    vector mpLocalCoords = msh[l].points().cellCoordinates(mp, mpIndex, true);

                    if (mpLocalCoords == -vector::one)
                    {
                        FatalError
                            << "Interpolation error at level " << l
                            << " and direction " << d
                            << ". Mirror point: " << mp
                            << " and colocated cell index: "
                            << mpIndex
                            << endl;
                        FatalError.exit();
                    }

                    // Index of <MeshType> left-bottom-aft cell w.r.t. mp
                    labelVector mpLBA = mpIndex;

                    if
                    (
                        (word(MeshType::typeName) == "colocated" ? true : (d != 0))
                        && (mpLocalCoords.x() < 0.5)
                    )
                    {
                        mpLBA.x() -= 1;
                    }
                    if
                    (
                        (word(MeshType::typeName) == "colocated" ? true : (d != 1))
                        && (mpLocalCoords.y() < 0.5)
                    )
                    {
                        mpLBA.y() -= 1;
                    }
                    if
                    (
                        (word(MeshType::typeName) == "colocated" ? true : (d != 2))
                        && (mpLocalCoords.z() < 0.5)
                    )
                    {
                        mpLBA.z() -= 1;
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
                    const vector v(interpolationWeights(mp,interpPoints,true));

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
                            << "Interpolation error at level " << l
                            << " and direction " << d
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

                x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k)
                    + omega * (H*2.0*boundaryValues_[d] - mpValue);
            }
        }
    }
}

}

}

}
