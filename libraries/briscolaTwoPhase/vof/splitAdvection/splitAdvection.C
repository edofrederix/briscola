#include "splitAdvection.H"

#include "addToRunTimeSelectionTable.H"
#include "rectilinearMesh.H"
#include "truncatedHex.H"
#include "truncatedPiped.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(splitAdvection, 0);
addToRunTimeSelectionTable(vof, splitAdvection, dictionary);

void splitAdvection::updateFlux
(
    const colocatedLowerFaceScalarField& phi,
    const label d
)
{
    const colocatedVectorField& n = normal_;

    const colocatedScalarField& cv =
        fvMsh_.template metrics<colocated>().cellVolumes();
    const colocatedVertexVectorField& v =
        fvMsh_.template metrics<colocated>().vertexCenters();

    const scalar dt = fvMsh_.time().deltaT().value();

    forAllFaces(phi, fd, i, j, k)
    if (d == fd)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk+units[d]);

        // For donating cells only

        if
        (
            (phi(ijk)[d] > 0 && alpha_(ijk) > vof::threshold)
         || (phi(ijk)[d] < 0 && alpha_(nei) > vof::threshold)
        )
        {
            labelVector don = phi(ijk)[d] > 0 ? ijk : nei;

            const scalar frac = phi(ijk)[d]*dt/cv(don);

            if (alpha_(don) > 1-vof::threshold)
            {
                flux_(ijk)[d] = cv(don)*frac/dt;
            }
            else if (Foam::mag(n(don)))
            {
                // Truncate the hexahedron in face-normal direction. This is not
                // necessarily volume conserving on general hexahedrons and
                // should be improved.

                const scalar C = lve_(alpha_(don),v(don),cv(don),n(don));

                const label l = phi(ijk)[d] > 0 ? d*2 : d*2+1;
                const label u = phi(ijk)[d] < 0 ? d*2 : d*2+1;

                vertexVector vertices(v(don));

                if (rectilinear_)
                {
                    for (int iv = 0; iv < 4; iv++)
                    {
                        vertices[vertexNumsInFace[u][iv]] =
                            vertices[vertexNumsInFace[l][iv]]
                          + frac
                          * (
                                vertices[vertexNumsInFace[u][iv]]
                              - vertices[vertexNumsInFace[l][iv]]
                            );
                    }
                }
                else
                {
                    scalar cutC = lve_.fluxVolumeLVE
                    (
                        v(don),
                        cv(don),
                        cv(don)*frac,
                        l,
                        u
                    );

                    for (int iv = 0; iv < 4; iv++)
                    {
                        vertices[vertexNumsInFace[u][iv]] =
                            vertices[vertexNumsInFace[u][iv]]
                          + cutC
                          * (
                                vertices[vertexNumsInFace[l][iv]]
                              - vertices[vertexNumsInFace[u][iv]]
                            );
                    }

                }

                scalar fluxVolume =
                    rectilinear_
                    ? truncatedPiped(vertices,n(don),C).volume()
                    : truncatedHex(vertices,n(don),C).volume();

                flux_(ijk)[d] = fluxVolume/dt;
            }
            else
            {
                flux_(ijk)[d] = 0.0;
            }
        }
    }
}

splitAdvection::splitAdvection
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    colocatedScalarField& alpha
)
:
    geometricVof(fvMsh, dict, normal, alpha),
    flux_
    (
        IOobject::groupName("flux", alpha.name()),
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    )
{}

splitAdvection::splitAdvection(const splitAdvection& vf)
:
    geometricVof(vf),
    flux_
    (
        IOobject::groupName("flux", vf.alpha().name()),
        vf.fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    )
{}

splitAdvection::~splitAdvection()
{}

void splitAdvection::solve(const colocatedLowerFaceScalarField& phi)
{
    alpha_.setOldTime();

    colocatedLabelField alphac
    (
        "alphac",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const colocatedScalarField& cv =
        fvMsh_.template metrics<colocated>().cellVolumes();

    // Central volume fraction value of Weymouth & Yue (2010)

    forAllCells(alpha_, i, j, k)
        alphac(i,j,k) = alpha_(i,j,k) > 0.5;

    // Solve the advection equation in a split way. Rotate directions.

    this->resetFlux();

    const scalar dt = fvMsh_.time().deltaT().value();
    const label ti = fvMsh_.time().timeIndex();

    for (int d = 0; d < 3; d++)
    {
        const label dir = (ti+d)%3;

        this->updateFlux(phi,dir);

        forAllCells(alpha_, i, j, k)
        {
            labelVector ijk(i,j,k);
            labelVector nei(ijk+units[dir]);

            alpha_(ijk) +=
                dt/cv(ijk)
              * (
                    scalar(alphac(ijk))*(phi(ijk)[dir] + phi(nei)[dir])
                  - flux_(ijk)[dir]
                  - flux_(nei)[dir]
                );
        }

        alpha_[0].correctBoundaryConditions();

        // Update the normal after the alpha update, so that it is consistent
        // with alpha and can be reused by other parts of the code

        normal_.correct();
    }
}

}

}

}
