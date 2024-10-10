// Velocity

colocatedVectorField U
(
    "U",
    fvMsh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE,
    true,
    true
);

U = Zero;

// Flux

colocatedLowerFaceScalarField phi
(
    "phi",
    fvMsh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE,
    true
);

phi = Zero;

// Curvature scheme

autoPtr<curvatureScheme> curvatureSchemePtr
(
    curvatureScheme::New
    (
        fvMsh,
        vofDict.subDict("curvatureScheme"),
        normal,
        alpha
    ).ptr()
);

curvatureScheme& curvature = curvatureSchemePtr();

// Write normal and kappa

normal.writeOpt() = IOobject::AUTO_WRITE;
curvature.writeOpt() = IOobject::AUTO_WRITE;
