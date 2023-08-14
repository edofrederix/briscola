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

const scalar vof::threshold_ = 1e-6;

void vof::updateFlux
(
    const colocatedFaceScalarField& phi,
    const label d
)
{
    const colocatedVectorField normal(this->normal()());

    const colocatedScalarDirection& a = alpha_[0][0];
    const colocatedVectorDirection& n = normal[0][0];
    const colocatedFaceScalarDirection& p = phi[0][0];
    const colocatedScalarDirection& V =
        fvMsh_.template metrics<colocated>().cellVolumes()[0][0];
    const colocatedVertexVectorDirection& v =
        fvMsh_.template metrics<colocated>().vertexCenters()[0][0];

    colocatedFaceScalarDirection& f = flux_[0][0];

    const scalar dt = fvMsh_.time().deltaT().value();

    forAllCells(f, i, j, k)
    {
        const labelVector ijk(i,j,k);

        const labelVector ijku(ijk + units[d]);
        const labelVector ijkl(ijk - units[d]);

        // For donating cells only

        if ((p(ijk)[d*2] > 0 || p(ijk)[d*2+1] > 0) && a(ijk) > threshold_)
        {
            const scalar magn = Foam::mag(n(ijk));

            const scalar C = magn > 0 ? lve_(ijk,n(ijk)) : 0;

            for (int fi = 0; fi < 2; fi++)
            if (p(ijk)[d*2+fi] > 0)
            {
                const label l = d*2 + fi;
                const label u = d*2 + (fi == 0);

                const scalar frac = p(ijk)[l]*dt/V(ijk);

                scalar flux = 0;

                if (a(ijk) > 1-threshold_)
                {
                    flux = V(ijk)*frac/dt;
                }
                else if (magn)
                {
                    // Truncate the hexahedron in face-normal direction.
                    // This is not necessarily volume conserving on general
                    // hexahedrons and should be improved.

                    vertexVector vertices(v(ijk));

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

                    scalar fluxVolume =
                        rectilinear_
                      ? truncatedPiped(vertices,n(ijk),C).volume()
                      : truncatedHex(vertices,n(ijk),C).volume();

                    flux = fluxVolume/dt;
                }

                if (flux != 0)
                {
                    f(ijk)[l] = flux;

                    // Also set the flux on the receiving end

                    if (fi == 0)
                    {
                        f(ijkl)[u] = -flux;
                    }
                    else
                    {
                        f(ijku)[u] = -flux;
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

            const labelVector S(f.boundaryStart(bo));
            const labelVector E(f.boundaryEnd(bo));

            labelVector ijk;

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                if (f(ijk+bo)[u] > 0)
                    f(ijk)[l] = -f(ijk+bo)[u];
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
        true
    ),
    flux_
    (
        IOobject::groupName("flux", alpha_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    rectilinear_(fvMsh.msh()[0].rectilinear() == unitXYZ),
    lve_(*this),
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

    const colocatedFaceScalarDirection& p = phi[0][0];
    const colocatedScalarDirection& V =
        fvMsh_.template metrics<colocated>().cellVolumes()[0][0];

    colocatedScalarDirection& a = alpha_[0][0];
    colocatedLabelDirection& ac = alphac[0][0];

    // Central volume fraction value of Weymouth & Yue (2010)

    forAllCells(ac, i, j, k)
        ac(i,j,k) = a(i,j,k) > 0.5;

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

        const colocatedFaceScalarDirection& f = flux_[0][0];

        forAllCells(a, i, j, k)
        {
            a(i,j,k) +=
                dt/V(i,j,k)
              * (
                    scalar(ac(i,j,k))*(p(i,j,k)[l] + p(i,j,k)[u])
                  - f(i,j,k)[l]
                  - f(i,j,k)[u]
                );
        }

        alpha_.correctBoundaryConditions();
    }
}

}

}

}
