// Velocity

staggeredScalarField U
(
    "U",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

U = Zero;
U.correctBoundaryConditions();

// Reconstruct colocated velocity field

colocatedVectorField Uc
(
    "Uc",
    fvMsh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE,
    true
);

Uc = Zero;

Switch colocatedReconstruction
(
    solverDict.lookupOrDefault
    (
        "colocatedReconstruction",
        false
    )
);

// Pressure

colocatedScalarField p
(
    "p",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

p = Zero;
p.correctBoundaryConditions();

// Extrapolated pressure (Dodd & Ferrante)

colocatedScalarField q
(
    "q",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true,
    false
);

q = Zero;

// Flux

staggeredLowerFaceScalarField phi
(
    "phi",
    fvMsh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE,
    true
);

phi = ex::faceFlux(U);

staggeredScalarField exSource
(
    "source",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

staggeredScalarField imSourceCoeff
(
    "sourceCoeff",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

exSource = scalarList(solverDict.lookup("source"));
imSourceCoeff = scalarList(solverDict.lookup("sourceCoeff"));

// Solvers

linearSystem<stencil,scalar,staggered> USys(U);
linearSystem<stencil,scalar,colocated> pSys(p);

autoPtr<solver<stencil,scalar,staggered>> USolve
(
    solver<stencil,scalar,staggered>::New
    (
        USys.x().name(),
        fvMsh
    )
);

autoPtr<PoissonSolver<stencil,scalar,colocated>> Poisson
(
    PoissonSolver<stencil,scalar,colocated>::New
    (
        pSys.x().name(),
        fvMsh
    )
);

// Specific volumes

staggeredScalarField v
(
    "v",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

v.setRestrictionScheme("harmonicVolumeWeighted");

colocatedScalarField vc
(
    "v",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

vc.setRestrictionScheme("harmonicVolumeWeighted");

colocatedLowerFaceScalarField vcf
(
    "v",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

v = icoTwoPhase.v<staggered>();
vc = icoTwoPhase.v<colocated>();
vcf = ex::interp(vc);

// Adams-Bashforth source

staggeredScalarField H
(
    "H",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true,
    false
);

H = ex::div(phi,U)
  - stagDotProduct(ex::grad(mu),ex::grad(U))*v;

// It is assumed that the minimum mass density is constant

const scalar minRho = gMin(icoTwoPhase.minRho<colocated>())[0];
const scalar maxv = 1.0/minRho;
