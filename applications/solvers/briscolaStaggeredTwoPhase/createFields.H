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

colocatedScalarField extrapolatedP
(
    "extrapolatedP",
    fvMsh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true,
    false
);

extrapolatedP = Zero;
extrapolatedP.correctBoundaryConditions();

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

// Viscous term

linearSystem<stencil,scalar,staggered> D(U);

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

H = ex::div(phi,U) - ((ex::grad(mu) & ex::grad(U)) / rho);