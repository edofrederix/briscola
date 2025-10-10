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

    if (!fvMsh.db().foundObject<colocatedFaceVectorField>("faceDotGradDelta0"))
    {
        // Cache fields

        tmp<colocatedFaceVectorField> tDelta0
        (
            new colocatedFaceVectorField
            (
                "faceDotGradDelta0",
                fvMsh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            )
        );

        tmp<colocatedFaceVectorField> tDelta1
        (
            new colocatedFaceVectorField
            (
                "faceDotGradDelta1",
                fvMsh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            )
        );

        tmp<colocatedFaceVectorField> tDelta2
        (
            new colocatedFaceVectorField
            (
                "faceDotGradDelta2",
                fvMsh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            )
        );

        colocatedFaceVectorDirection& delta0 = tDelta0.ref()[0][0];
        colocatedFaceVectorDirection& delta1 = tDelta1.ref()[0][0];
        colocatedFaceVectorDirection& delta2 = tDelta2.ref()[0][0];

        const meshDirection<faceVector,colocated>& fn =
            fvMsh.template metrics<colocated>().faceNormals()[0][0];

        const meshDirection<faceScalar,colocated>& delta =
            fvMsh.template metrics<colocated>().faceDeltas()[0][0];

        const meshDirection<faceVector,colocated>& fc =
            fvMsh.template metrics<colocated>().faceCenters()[0][0];

        // Normal direction coefficient. The minus sign is because the face
        // normal points in negative direction on the lower face

        forAllFacesInDirection(delta0, fd, i, j, k)
            delta0(i,j,k)[fd*2] = -delta(i,j,k)[fd*2]*fn(i,j,k)[fd*2];

        // Tangential direction coefficients

        forAllFacesInDirection(delta1, fd, i, j, k)
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

            delta1(i,j,k)[fd*2] = dist1/Foam::magSqr(dist1);
            delta2(i,j,k)[fd*2] = dist2/Foam::magSqr(dist2);
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

    tmp<meshField<FaceSpace<Type>,colocated>> tGrad
    (
        new meshField<FaceSpace<Type>,colocated>
        (
            "faceDotGrad("+field.name()+")",
            fvMsh
        )
    );

    meshDirection<FaceSpace<Type>,colocated>& grad = tGrad.ref()[0][0];

    // Cached coefficients

    cache(fvMsh);

    const colocatedFaceVectorDirection& delta0 =
        fvMsh.db()
       .lookupObject<colocatedFaceVectorField>("faceDotGradDelta0")[0][0];

    const colocatedFaceVectorDirection& delta1 =
        fvMsh.db()
       .lookupObject<colocatedFaceVectorField>("faceDotGradDelta1")[0][0];

    const colocatedFaceVectorDirection& delta2 =
        fvMsh.db()
       .lookupObject<colocatedFaceVectorField>("faceDotGradDelta2")[0][0];

    const meshDirection<faceVector,colocated>& fan =
        fvMsh.template metrics<colocated>().faceAreaNormals()[0][0];

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

        grad(ijk)[fd*2] =  delta0(ijk)[fd*2]*(diff0 & fan(ijk)[fd*2]);
        grad(ijk)[fd*2] += delta1(ijk)[fd*2]*(diff1 & fan(ijk)[fd*2]);
        grad(ijk)[fd*2] += delta2(ijk)[fd*2]*(diff2 & fan(ijk)[fd*2]);

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

    tmp<meshField<FaceSpace<Type>,staggered>> tGrad
    (
        new meshField<FaceSpace<Type>,staggered>
        (
            "faceDotGrad("+field.name()+")",
            fvMsh
        )
    );

    meshLevel<FaceSpace<Type>,staggered>& grad = tGrad.ref()[0];

    const meshLevel<faceScalar,staggered>& delta =
        fvMsh.template metrics<staggered>().faceDeltas()[0];

    const meshLevel<faceScalar,staggered>& fa =
        fvMsh.template metrics<staggered>().faceAreas()[0];

    forAllFaces(grad, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        // The gradient in d-direction of fd-direction staggered components

        grad(d,ijk)[fd*2] =
            delta(fd,ijk)[d*2]*fa(fd,ijk)[d*2]
          * (field(fd,lowerNeighbor(ijk,d)) - field(fd,ijk));

        grad(d,nei)[fd*2+1] = -grad(d,ijk)[fd*2];
    }

    return tGrad;
}

}

}

}
