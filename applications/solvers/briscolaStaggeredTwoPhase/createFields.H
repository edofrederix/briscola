// Velocity

staggeredScalarField U
(
    "U",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true
);

U = Zero;
U.correctBoundaryConditions();

// Reconstructed colocated velocity

colocatedVectorField Uc
(
    "Uc",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

Uc = Zero;

// Pressure

colocatedScalarField p
(
    "p",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true
);

p = Zero;
p.correctBoundaryConditions();

// Flux

staggeredFaceScalarField phi
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

linearSystem<stencil,scalar,staggered> USys("NS(U)", U, true);

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
        p.name(),
        fvMsh
    )
);

Poisson->enableFlux();

// Runge-Kutta stage sources

PtrList<staggeredScalarField> stageSourcesA
(
    rk.stageList<staggeredScalarField>("A")
);

PtrList<staggeredScalarField> stageSourcesB
(
    rk.stageList<staggeredScalarField>("B")
);

linearSystem<stencil,scalar,staggered> USysA("A(U)", U);
linearSystem<stencil,scalar,staggered> USysB("B(U)", U);

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

colocatedFaceScalarField vcf
(
    "vcf",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
);

vcf.setRestrictionScheme("faceAreaWeighted");
