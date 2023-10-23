#include "vof.H"
#include "rectilinearMesh.H"
#include "truncatedHex.H"
#include "truncatedPiped.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(vof, 0);

const scalar vof::threshold = 1e-12;

void vof::updateFlux
(
    const colocatedFaceScalarField& phi,
    const label d
)
{
    const colocatedVectorField n(this->normal()());

    const colocatedScalarField& cv =
        fvMsh_.template metrics<colocated>().cellVolumes();
    const colocatedVertexVectorField& v =
        fvMsh_.template metrics<colocated>().vertexCenters();

    const scalar dt = fvMsh_.time().deltaT().value();

    forAllCells(alpha_, i, j, k)
    {
        const labelVector ijk(i,j,k);

        const labelVector ijku(ijk + units[d]);
        const labelVector ijkl(ijk - units[d]);

        // For donating cells only

        if
        (
            (phi(ijk)[d*2] > 0 || phi(ijk)[d*2+1] > 0)
         && alpha_(ijk) > vof::threshold
        )
        {
            const scalar magn = Foam::mag(n(ijk));

            const scalar C =
                magn > 0 ? lve_(alpha_(ijk),v(ijk),cv(ijk),n(ijk)) : 0;

            for (int fi = 0; fi < 2; fi++)
            if (phi(ijk)[d*2+fi] > 0)
            {
                const label l = d*2 + fi;
                const label u = d*2 + (fi == 0);

                const scalar frac = phi(ijk)[l]*dt/cv(ijk);

                scalar flux = 0;

                if (alpha_(ijk) > 1-vof::threshold)
                {
                    flux = cv(ijk)*frac/dt;
                }
                else if (magn)
                {
                    // Truncate the hexahedron in face-normal direction. This is
                    // not necessarily volume conserving on general hexahedrons
                    // and should be improved.

                    vertexVector vertices(v(ijk));

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
                            v(ijk),
                            cv(ijk),
                            frac * cv(ijk),
                            u,
                            l
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
                      ? truncatedPiped(vertices,n(ijk),C).volume()
                      : truncatedHex(vertices,n(ijk),C).volume();

                    flux = fluxVolume/dt;
                }

                if (flux != 0)
                {
                    flux_(ijk)[l] = flux;

                    // Also set the flux on the receiving end

                    if (fi == 0)
                    {
                        flux_(ijkl)[u] = -flux;
                    }
                    else
                    {
                        flux_(ijku)[u] = -flux;
                    }
                }
            }
        }
    }

    // Add the ghost boundary fluxes to receiving boundary cells

    flux_.correctCommBoundaryConditions();

    forAll(flux_.boundaryConditions(), bci)
    {
        const auto& bc = flux_.boundaryConditions()[bci];

        if
        (
            bc.boundaryOffsetDegree() == 1
         && (
                bc.baseType() == PARALLELBC
             || bc.baseType() == PERIODICBC
            )
        )
        {
            const labelVector bo = bc.boundaryOffset();

            const label l = faceNumber( bo);
            const label u = faceNumber(-bo);

            const labelVector S(fvMsh_.S<colocated>(bo));
            const labelVector E(fvMsh_.E<colocated>(bo));

            labelVector ijk;

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                if (flux_(ijk+bo)[u] > 0)
                    flux_(ijk)[l] = -flux_(ijk+bo)[u];
            }
        }
    }
}

vof::vof(const IOdictionary& dict, const fvMesh& fvMsh)
:
    regIOobject(dict, true),
    fvMsh_(fvMsh),
    alpha_
    (
        "alpha",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        true,
        true,
        false
    ),
    flux_
    (
        IOobject::groupName("flux", alpha_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        false
    ),
    rectilinear_(fvMsh.msh()[0].rectilinear() == unitXYZ),
    lve_(rectilinear_),
    normalSchemePtr_
    (
        normalScheme::New
        (
            *this,
            dict.subDict("normalScheme")
        ).ptr()
    )
{}

vof::~vof()
{}

void vof::solve(const colocatedFaceScalarField& phi)
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
        const label l = dir*2;
        const label u = dir*2+1;

        this->updateFlux(phi,dir);

        forAllCells(alpha_, i, j, k)
        {
            alpha_(i,j,k) +=
                dt/cv(i,j,k)
              * (
                    scalar(alphac(i,j,k))*(phi(i,j,k)[l] + phi(i,j,k)[u])
                  - flux_(i,j,k)[l]
                  - flux_(i,j,k)[u]
                );
        }

        alpha_.correctBoundaryConditions();
    }
}

}

}

}
