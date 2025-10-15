#include "linearFaceDotGradientScheme.H"
#include "linearInterpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
void linearFaceDotGradientScheme<Type,colocated>::cache
(
    const fvMesh& fvMsh
) const
{
    // Create coefficients and store

    if (!fvMsh.db().foundObject<staggeredVectorField>("faceDotGradDelta0"))
    {
        // Cache fields

        tmp<staggeredVectorField> tDelta0 =
            staggeredVectorField::New
            (
                "faceDotGradDelta0",
                fvMsh
            );

        tmp<staggeredVectorField> tDelta1 =
            staggeredVectorField::New
            (
                "faceDotGradDelta1",
                fvMsh
            );

        tmp<staggeredVectorField> tDelta2 =
            staggeredVectorField::New
            (
                "faceDotGradDelta2",
                fvMsh
            );

        staggeredVectorLevel& delta0 = tDelta0.ref()[0];
        staggeredVectorLevel& delta1 = tDelta1.ref()[0];
        staggeredVectorLevel& delta2 = tDelta2.ref()[0];

        const FastPtrList<colocatedVectorField>& fn =
            fvMsh.template metrics<colocated>().soa().faceNormals();

        const FastPtrList<colocatedScalarField>& delta =
            fvMsh.template metrics<colocated>().soa().faceDeltas();

        const FastPtrList<colocatedVectorField>& fc =
            fvMsh.template metrics<colocated>().soa().faceCenters();

        // Normal direction coefficient. The minus sign is because the face
        // normal points in negative direction on the lower face

        forAllFacesInDirection(delta[0], fd, i, j, k)
            delta0(fd,i,j,k) = -delta[fd](i,j,k)*fn[fd](i,j,k);

        // Tangential direction coefficients

        forAllFacesInDirection(delta[0], fd, i, j, k)
        {
            const labelVector ijk(i,j,k);
            const labelVector nei(lowerNeighbor(ijk,fd));

            // Indices of the other two direction

            const label a = fd == 0 ? 1 : 0;
            const label b = fd == 2 ? 1 : 2;

            // Tangential distances across the face

            const vector dist1 =
                fc[fd](upperNeighbor(ijk,a))
              - fc[fd](lowerNeighbor(ijk,a));

            const vector dist2 =
                fc[fd](upperNeighbor(ijk,b))
              - fc[fd](lowerNeighbor(ijk,b));

            delta1(fd,i,j,k) = dist1/Foam::magSqr(dist1);
            delta2(fd,i,j,k) = dist2/Foam::magSqr(dist2);
        }

        // Store

        fvMsh.db().objectRegistry::store(tDelta0.ptr());
        fvMsh.db().objectRegistry::store(tDelta1.ptr());
        fvMsh.db().objectRegistry::store(tDelta2.ptr());
    }
}

template<class Type>
tmp<meshField<FaceSpace<Type>,colocated>>
linearFaceDotGradientScheme<Type,colocated>::faceDotGrad
(
    const meshField<Type,colocated>& field
)
{
    const fvMesh& fvMsh = field.fvMsh();

    tmp<meshField<FaceSpace<Type>,colocated>> tGrad =
        meshField<FaceSpace<Type>,colocated>::New
        (
            "faceDotGrad("+field.name()+")",
            fvMsh
        );

    meshDirection<FaceSpace<Type>,colocated>& grad = tGrad.ref()[0][0];

    // Cached coefficients

    cache(fvMsh);

    const staggeredVectorLevel& delta0 =
        fvMsh.db().lookupObject<staggeredVectorField>("faceDotGradDelta0")[0];

    const staggeredVectorLevel& delta1 =
        fvMsh.db().lookupObject<staggeredVectorField>("faceDotGradDelta1")[0];

    const staggeredVectorLevel& delta2 =
        fvMsh.db().lookupObject<staggeredVectorField>("faceDotGradDelta2")[0];

    const FastPtrList<colocatedVectorField>& fan =
        fvMsh.template metrics<colocated>().soa().faceAreaNormals();

    // Compute face-dot gradient

    forAllFacesInDirection(grad, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(ijk,fd));

        // Indices of the other two direction

        const label a = fd == 0 ? 1 : 0;
        const label b = fd == 2 ? 1 : 2;

        // Normal difference from cell centered values

        const Type diff0 = field(ijk) - field(nei);

        // Tangential differences from values interpolated to the face (could be
        // improved with face weights)

        const Type diff1 =
            0.5*(field(upperNeighbor(ijk,a)) + field(upperNeighbor(nei,a)))
          - 0.5*(field(lowerNeighbor(ijk,a)) + field(lowerNeighbor(nei,a)));

        const Type diff2 =
            0.5*(field(upperNeighbor(ijk,b)) + field(upperNeighbor(nei,b)))
          - 0.5*(field(lowerNeighbor(ijk,b)) + field(lowerNeighbor(nei,b)));

        // Gradient contributions

        grad(ijk)[fd*2] =  delta0(fd,ijk)*(diff0 & fan[fd](ijk));
        grad(ijk)[fd*2] += delta1(fd,ijk)*(diff1 & fan[fd](ijk));
        grad(ijk)[fd*2] += delta2(fd,ijk)*(diff2 & fan[fd](ijk));

        // Copy to neighbor

        grad(nei)[fd*2+1] = -grad(ijk)[fd*2];
    }

    return tGrad;
}

// Staggered

template<class Type>
tmp<meshField<FaceSpace<Type>,staggered>>
linearFaceDotGradientScheme<Type,staggered>::faceDotGrad
(
    const meshField<Type,staggered>& field
)
{
    const fvMesh& fvMsh = field.fvMsh();

    tmp<meshField<FaceSpace<Type>,staggered>> tGrad =
        meshField<FaceSpace<Type>,staggered>::New
        (
            "faceDotGrad("+field.name()+")",
            fvMsh
        );

    meshLevel<FaceSpace<Type>,staggered>& grad = tGrad.ref()[0];

    const FastPtrList<staggeredScalarField>& delta =
        fvMsh.template metrics<staggered>().soa().faceDeltas();

    const FastPtrList<staggeredScalarField>& fa =
        fvMsh.template metrics<staggered>().soa().faceAreas();

    forAllFaces(grad, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        // The gradient in d-direction of fd-direction staggered components

        grad(d,ijk)[fd*2] =
            delta[d](fd,ijk)*fa[d](fd,ijk)
          * (field(fd,lowerNeighbor(ijk,d)) - field(fd,ijk));

        grad(d,nei)[fd*2+1] = -grad(d,ijk)[fd*2];
    }

    return tGrad;
}

}

}

}
