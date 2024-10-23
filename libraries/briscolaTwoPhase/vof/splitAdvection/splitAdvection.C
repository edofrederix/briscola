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
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNei(i,j,k,d));

        // Select the donating side

        const scalar flux = phi(ijk)[d];
        const labelVector don(flux > 0 ? ijk : nei);

        if (Foam::mag(flux) > 0 && alpha_(don) > vof::threshold)
        {
            const scalar frac = Foam::mag(flux)*dt/cv(don);

            scalar fluxAlpha = 0;

            if (alpha_(don) >= 1 - vof::threshold)
            {
                fluxAlpha = cv(don)*frac/dt;
            }
            else if (Foam::mag(n(don)))
            {
                const scalar C = lve_(alpha_(don),v(don),cv(don),n(don));

                const label l = flux > 0 ? d*2 : d*2+1;
                const label u = flux < 0 ? d*2 : d*2+1;

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
                            vertices[vertexNumsInFace[l][iv]]
                          + cutC
                          * (
                                vertices[vertexNumsInFace[u][iv]]
                              - vertices[vertexNumsInFace[l][iv]]
                            );
                    }

                }

                scalar fluxVolume =
                    rectilinear_
                    ? truncatedPiped(vertices,n(don),C).volume()
                    : truncatedHex(vertices,n(don),C).volume();

                fluxAlpha = fluxVolume/dt;
            }
            else
            {
                // We have an interfacial cell without normal

                #ifdef FULLDEBUG

                WarningInFunction
                    << "Donating cell has interface but no normal" << nl;

                #endif
            }

            flux_(ijk)[d] = Foam::sign(flux)*fluxAlpha;
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

    const faceLabel boundaryType(this->fvMsh_.msh().faceBoundaryType());

    for (int d = 0; d < 3; d++)
    {
        const label dir = (ti + d) % 3;

        // Skip empty directions

        if (boundaryType[dir*2] != emptyBoundary::typeNumber)
        {
            this->updateFlux(phi, dir);

            forAllCells(alpha_, i, j, k)
            {
                labelVector ijk(i,j,k);
                labelVector nei(upperNei(i,j,k,dir));

                alpha_(ijk) +=
                    dt/cv(ijk)
                  * (
                        scalar(alphac(ijk))*(phi(ijk)[dir] - phi(nei)[dir])
                      - flux_(ijk)[dir]
                      + flux_(nei)[dir]
                    );

                // Remove tiny errors in alpha

                alpha_(ijk) = Foam::min(Foam::max(alpha_(ijk), 0.0), 1.0);
            }

            alpha_.correctBoundaryConditions();

            // Update the normal after the alpha update, so that it is
            // consistent with alpha and can be reused by other parts of the
            // code

            normal_.correct();
        }
    }
}

}

}

}
