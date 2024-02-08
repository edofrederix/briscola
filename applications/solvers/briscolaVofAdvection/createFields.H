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

autoPtr<curvatureScheme> curvatureSchemePtr;

if (curvature)
    curvatureSchemePtr =
        curvatureScheme::New
        (
            fvMsh,
            vofDict.subDict("curvatureScheme"),
            normal,
            alpha
        ).ptr();

// Write normal and kappa

normal.writeOpt() = IOobject::AUTO_WRITE;

if (curvature)
    curvatureSchemePtr->writeOpt() = IOobject::AUTO_WRITE;
