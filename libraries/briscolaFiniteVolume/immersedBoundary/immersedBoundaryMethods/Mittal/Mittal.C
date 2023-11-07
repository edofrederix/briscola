#include "Mittal.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
Mittal<Type,MeshType>::Mittal
(
    dictionary& dict,
    const fvMesh& fvMsh
)
:
    immersedBoundaryMethod<Type,MeshType>(dict,fvMsh,true)
{

}

// Destructor

template<class Type, class MeshType>
Mittal<Type,MeshType>::~Mittal()
{}

template<class Type, class MeshType>
void Mittal<Type,MeshType>::correctJacobiPoints
(
    meshLevel<Type,MeshType>& x
)
{
    scalar omega = this->omega_;

    // Mesh
    const fvMesh& fvMsh = x.fvMsh();
    const mesh& msh = fvMsh.msh();

    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh.metrics<MeshType>().cellCenters();

    scalar tol = 1e-5;

    label l = x.levelNum();

    forAllDirections(x,d,i,j,k)
    {
        if (this->ghostMask_(l,d,i,j,k) == 1)
        {
            // mirror point - make this a field in IB
            vector mp = this->mirrorPoint(CC(l,d,i,j,k));

            // Fix situations where the mirror point is just outside of the mesh
            // bounding box due to rounding errors
            if
            (
                (mp.x() <= msh[l].boundingBox().left() + tol)
                && (mp.x() >= msh[l].boundingBox().left() - tol)
            )
            {
                mp.x() += tol;
            }
            if
            (
                (mp.x() >= msh[l].boundingBox().right() - tol)
                && (mp.x() <= msh[l].boundingBox().right() + tol)
            )
            {
                mp.x() -= tol;
            }
            if
            (
                (mp.y() <= msh[l].boundingBox().bottom() + tol)
                && (mp.y() >= msh[l].boundingBox().bottom() - tol)
            )
            {
                mp.y() += tol;
            }
            if
            (
                (mp.y() >= msh[l].boundingBox().top() - tol)
                && (mp.y() <= msh[l].boundingBox().top() + tol)
            )
            {
                mp.y() -= tol;
            }
            if
            (
                (mp.z() <= msh[l].boundingBox().aft() + tol)
                && (mp.z() >= msh[l].boundingBox().aft() - tol)
            )
            {
                mp.z() += tol;
            }
            if
            (
                (mp.z() >= msh[l].boundingBox().fore() - tol)
                && (mp.z() <= msh[l].boundingBox().fore() + tol)
            )
            {
                mp.z() -= tol;
            }

            // Colocated cell index of mp
            labelVector mpIndex = msh.findCell(mp, l);

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

            // Index of staggered left-bottom-aft cell w.r.t. mp
            labelVector mpLBA = mpIndex;
            if
            (
                (d != 0)
                && (mpLocalCoords.x() < 0.5)
            )
            {
                mpLBA.x() -= 1;
            }
            if
            (
                (d != 1)
                && (mpLocalCoords.y() < 0.5)
            )
            {
                mpLBA.y() -= 1;
            }
            if
            (
                (d != 2)
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

            Type mpValue = Zero;

            if (v != -vector::one)
            {
                mpValue =
                      weights.lba()*x[d](mpLBA)
                    + weights.rba()*x[d](mpLBA+unitX)
                    + weights.lta()*x[d](mpLBA+unitY)
                    + weights.rta()*x[d](mpLBA+unitX+unitY)
                    + weights.lbf()*x[d](mpLBA+unitZ)
                    + weights.rbf()*x[d](mpLBA+unitX+unitZ)
                    + weights.ltf()*x[d](mpLBA+unitY+unitZ)
                    + weights.rtf()*x[d](mpLBA+unitX+unitY+unitZ);
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

            x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k) - omega * mpValue;
        }
    }
}

}

}

}
