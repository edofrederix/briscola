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

    if (!fvMsh.db().foundObject<colocatedVectorFaceField>("faceDotGradDelta0"))
    {
        // Cache fields

        tmp<colocatedVectorFaceField> tDelta0 =
            colocatedVectorFaceField::New
            (
                "faceDotGradDelta0",
                fvMsh
            );

        tmp<colocatedVectorFaceField> tDelta1 =
            colocatedVectorFaceField::New
            (
                "faceDotGradDelta1",
                fvMsh
            );

        tmp<colocatedVectorFaceField> tDelta2 =
            colocatedVectorFaceField::New
            (
                "faceDotGradDelta2",
                fvMsh
            );

        colocatedVectorFaceField& delta0 = tDelta0.ref();
        colocatedVectorFaceField& delta1 = tDelta1.ref();
        colocatedVectorFaceField& delta2 = tDelta2.ref();

        const colocatedVectorFaceField& fn =
            fvMsh.template metrics<colocated>().faceNormals();

        const colocatedScalarFaceField& delta =
            fvMsh.template metrics<colocated>().faceDeltas();

        // Use AoS storage for face centers because we need to access ghost
        // values. Because the result is cached this has no performance penalty.

        const colocatedFaceVectorField fc
        (
            fvMsh.template metrics<colocated>().aos().faceCenters()
        );

        // Normal direction coefficient. The minus sign is because the face
        // normal points in negative direction on the lower face

        forAllFaces(delta0, fd, i, j, k)
            delta0[fd](i,j,k) = -delta[fd](i,j,k)*fn[fd](i,j,k);

        // Tangential direction coefficients

        forAllFaces(delta1, fd, i, j, k)
        {
            const labelVector ijk(i,j,k);
            const labelVector nei(lowerNeighbor(ijk,fd));

            // Indices of the other two direction

            const label a = fd == 0 ? 1 : 0;
            const label b = fd == 2 ? 1 : 2;

            // Tangential distances across the face

            const vector dist1 =
                fc(upperNeighbor(ijk,a))[fd*2]
              - fc(lowerNeighbor(ijk,a))[fd*2];

            const vector dist2 =
                fc(upperNeighbor(ijk,b))[fd*2]
              - fc(lowerNeighbor(ijk,b))[fd*2];

            delta1[fd](i,j,k) = dist1/Foam::magSqr(dist1);
            delta2[fd](i,j,k) = dist2/Foam::magSqr(dist2);
        }

        // Store

        fvMsh.db().objectRegistry::store(tDelta0.ptr());
        fvMsh.db().objectRegistry::store(tDelta1.ptr());
        fvMsh.db().objectRegistry::store(tDelta2.ptr());
    }
}

template<class Type>
tmp<faceField<Type,colocated>>
linearFaceDotGradientScheme<Type,colocated>::faceDotGrad
(
    const meshField<Type,colocated>& field
)
{
    const fvMesh& fvMsh = field.fvMsh();

    tmp<faceField<Type,colocated>> tGrad =
        faceField<Type,colocated>::New
        (
            "faceDotGrad("+field.name()+")",
            fvMsh
        );

    faceField<Type,colocated>& grad = tGrad.ref();

    // Cached coefficients

    cache(fvMsh);

    const colocatedVectorFaceField& delta0 =
        fvMsh.db().lookupObject<colocatedVectorFaceField>("faceDotGradDelta0");

    const colocatedVectorFaceField& delta1 =
        fvMsh.db().lookupObject<colocatedVectorFaceField>("faceDotGradDelta1");

    const colocatedVectorFaceField& delta2 =
        fvMsh.db().lookupObject<colocatedVectorFaceField>("faceDotGradDelta2");

    const colocatedVectorFaceField& fan =
        fvMsh.template metrics<colocated>().faceAreaNormals();

    // Compute face-dot gradient

    forAllFaces(grad, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

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

        grad[fd](ijk) =
            delta0[fd](ijk)*(diff0 & fan[fd](ijk))
          + delta1[fd](ijk)*(diff1 & fan[fd](ijk))
          + delta2[fd](ijk)*(diff2 & fan[fd](ijk));
    }

    return tGrad;
}

// Staggered

template<class Type>
tmp<faceField<Type,staggered>>
linearFaceDotGradientScheme<Type,staggered>::faceDotGrad
(
    const meshField<Type,staggered>& field
)
{
    const fvMesh& fvMsh = field.fvMsh();

    tmp<faceField<Type,staggered>> tGrad =
        faceField<Type,staggered>::New
        (
            "faceDotGrad("+field.name()+")",
            fvMsh
        );

    faceField<Type,staggered>& grad = tGrad.ref();

    const staggeredScalarFaceField& delta =
        fvMsh.template metrics<staggered>().faceDeltas();

    const staggeredScalarFaceField& fa =
        fvMsh.template metrics<staggered>().faceAreas();

    forAllFaces(grad, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);

        // The gradient in d-direction of fd-direction staggered components

        grad[fd](d,ijk) =
            delta[d](fd,ijk)*fa[d](fd,ijk)
          * (field(fd,lowerNeighbor(i,j,k,d)) - field(fd,ijk));
    }

    return tGrad;
}

}

}

}
