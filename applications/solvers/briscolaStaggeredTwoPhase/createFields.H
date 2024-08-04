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

// Velocity equation residual

staggeredScalarField residual
(
    "residual",
    fvMsh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE,
    true
);

residual = Zero;

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

// Sources

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

exSource =
    scalarList
    (
        twoPhaseDict.lookupOrDefault
        (
            "source",
            scalarList(3, 0.0)
        )
    );

imSourceCoeff =
    scalarList
    (
        twoPhaseDict.lookupOrDefault
        (
            "sourceCoeff",
            scalarList(3, 0.0)
        )
    );

// Solvers

linearSystem<stencil,scalar,staggered> USys(U);
linearSystem<stencil,scalar,staggered> USysAux(U);
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

colocatedLowerFaceScalarField vcf
(
    "vcf",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

vcf.setRestrictionScheme("faceAreaWeighted");

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

H = Zero;

// It is assumed that the minimum mass density is constant

const scalar minRho = gMin(icoTwoPhase.minRho<colocated>())[0];
const scalar maxv = 1.0/minRho;

// Solver parameters

const Switch split = twoPhaseDict.lookupOrDefault<Switch>("split", false);
const Switch reduced = twoPhaseDict.lookupOrDefault<Switch>("reduced", false);
const label nCorr = twoPhaseDict.lookupOrDefault<label>("nCorrectors", 1);
