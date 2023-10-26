#include "incompressibleVof.H"

#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(incompressibleVof, 0);
addToRunTimeSelectionTable(twoPhaseModel, incompressibleVof, dictionary);

incompressibleVof::incompressibleVof
(
    const IOdictionary& dict,
    const fvMesh& fvMsh
)
:
    twoPhaseModel(dict, fvMsh),
    vfPtr_(vof::New(dict.subDict("vof"), fvMsh)),
    rho1_(readScalar(dict.lookup("rho1"))),
    rho2_(readScalar(dict.lookup("rho2"))),
    mu1_(readScalar(dict.lookup("mu1"))),
    mu2_(readScalar(dict.lookup("mu2"))),
    fluxName_(dict.lookupOrDefault<word>("flux", "phi")),
    velocityName_(dict.lookupOrDefault<word>("velocity", "U"))
{}

incompressibleVof::incompressibleVof(const incompressibleVof& tpm)
:
    twoPhaseModel(tpm),
    vfPtr_(tpm.vf().clone()),
    rho1_(tpm.rho1_),
    rho2_(tpm.rho2_),
    mu1_(tpm.mu1_),
    mu2_(tpm.mu2_),
    fluxName_(tpm.fluxName_),
    velocityName_(tpm.velocityName_)
{}

incompressibleVof::~incompressibleVof()
{}

void incompressibleVof::correct()
{
    // Update the volume fraction by solving the Vof advection equation

    if (fvMsh_.db().foundObject<colocatedFaceScalarField>(fluxName_))
    {
        const colocatedFaceScalarField& phi =
            fvMsh_.db().lookupObject<colocatedFaceScalarField>(fluxName_);

        vfPtr_->solve(phi);
    }
    else
    {
        const staggeredScalarField& U =
            fvMsh_.db().lookupObject<staggeredScalarField>(velocityName_);

        vfPtr_->solve(U);
    }

    // Correct the colocated mixture mass density and dynamic viscosity

    colocatedScalarField& alpha = vfPtr_->alpha();
    alpha.correctBoundaryConditions();

    const colocatedFaceScalarField alphaf
    (
        ex::interp(alpha)
    );

    rhoc_ = rho2_*alpha + rho1_*(1.0-alpha);
    muc_ = 1.0/(alphaf/mu2_ + (1.0-alphaf)/mu1_);

    // Correct the staggered mixture mass density and dynamic viscosity

    if (rhosPtr_.valid() && musPtr_.valid())
    {
        staggeredScalarField& rhos = rhosPtr_();
        staggeredFaceScalarField& mus = musPtr_();

        const staggeredScalarField alphas
        (
            stagInterp(alpha)
        );

        const staggeredFaceScalarField alphafs
        (
            stagFaceInterp(alpha)
        );

        rhos = rho2_*alphas + rho1_*(1.0-alphas);
        mus = 1.0/(alphafs/mu2_ + (1.0-alphafs)/mu1_);
    }
}

}

}

}
