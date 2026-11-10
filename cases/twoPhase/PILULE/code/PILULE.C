#include "PILULE.H"
#include "addToRunTimeSelectionTable.H"
#include "meshFields.H"
#include "List.H"
#include "PstreamReduceOps.H"
#include "Time.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(PILULE, 0);

addToRunTimeSelectionTable
(
    functionObject,
    PILULE,
    dictionary
);

PILULE::PILULE
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    fluidVol_(Zero),
    removedVol_(Zero),
    volLossIBM_
    (
        "volLossIBM",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true,
        false
    ),
    U1_
    (
        "U1",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false
    ),
    Usqr1_
    (
        "Usqr1",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false
    ),
    U2_
    (
        "U2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false
    ),
    Usqr2_
    (
        "Usqr2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false
    ),
    piluleDict_
    (
        IOobject
        (
            "piluleDict",
            runTime.time().constant(),
            runTime.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    oldError_(0),
    kP_(readScalar(dict.lookup("Kp"))),
    kI_(readScalar(dict.lookup("Ki"))),
    targetJ2_(readScalar(dict.lookup("Jliq"))),
    forcing_(0)
{
    volLossIBM_ = Zero;

    U1_ = Zero;
    U2_ = Zero;
    Usqr1_ = Zero;
    Usqr2_ = Zero;

    piluleDict_.add("alphaH", false);
    piluleDict_.store();

    piluleDict_.set("alphaH", false);

    read(dict);

    if (Pstream::master())
    {
        filePtr_.reset(new OFstream("data.txt"));
    }

    const colocatedScalarDirection& cv =
        fvMsh_.metrics<colocated>().cellVolumes()[0][0];

    const colocatedScalarDirection& ibmMask =
        fvMsh_.metrics<colocated>().mask()[0][0];

    fluidVol_ = gSum(cv*(1.0-ibmMask));

    piluleDict_.add("fluidVol", false);
    piluleDict_.set("fluidVol", fluidVol_);
}

PILULE::~PILULE()
{}

bool PILULE::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        // Initialize velocity

        const vector U0(0,0,0.05);

        if (runTime_.foundObject<staggeredScalarField>("U"))
        {
            // Staggered

            staggeredScalarField& U =
                runTime_.lookupObjectRef<staggeredScalarField>("U");

            const tensor base =
                U.fvMsh().msh().cast<rectilinearMesh>().base();

            forAllCells(U,d,i,j,k)
            {
                U(d,i,j,k) = (base & U0)[d];
            }

            U.correctBoundaryConditions();
        }
        else
        {
            // Colocated

            colocatedVectorField& U =
                runTime_.lookupObjectRef<colocatedVectorField>("U");

            forAllCells(U,i,j,k)
            {
                U(i,j,k) = U0;
            }

            U.correctBoundaryConditions();
        }

        // Initialize bubbles

        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = Zero;

        const colocatedVertexVectorField& v =
            alpha.fvMsh().metrics<colocated>().vertexCenters();

        vectorList bubbleCenters(0);

        const scalar bubbleRadius(0.0015);
        const scalar pipeRadius(0.019);
        const scalar pipeLength(0.16);
        const scalar cylinderRadius(0.006);
        // minimum spacing between bubbles or between bubble and wall
        const scalar minBubbleSpacing(1e-3);

        const label nBubbles(600);

        // Random bubble location seed
        std::srand(12345);

        label iter = 0;

        if (Pstream::master())
        {
            // Set bubbles randomly.

            while (bubbleCenters.size() < nBubbles)
            {
                scalar x = (2.0 * std::rand() / RAND_MAX - 1.0) * pipeRadius;
                scalar y = (2.0 * std::rand() / RAND_MAX - 1.0) * pipeRadius;
                scalar z = bubbleRadius
                    + std::rand() / (double)RAND_MAX
                    * (pipeLength - 2.0 * bubbleRadius);

                vector newBubble(x,y,z);

                // Check if bubble is inside IBM pipe walls
                if
                (
                    Foam::sqr(newBubble.x()) + Foam::sqr(newBubble.y()) <
                    Foam::sqr(pipeRadius-bubbleRadius-minBubbleSpacing)
                )
                {
                    bool overlap = false;

                    // Check whether there is overlap with cylinder
                    vector cylinderAxis(0,1,0);
                    vector AP = newBubble - vector(0,0,0.05);
                    vector proj =
                        ((AP & cylinderAxis) / magSqr(cylinderAxis))
                      * cylinderAxis;

                    scalar dist = mag(AP - proj);

                    if (dist < cylinderRadius + bubbleRadius + minBubbleSpacing)
                    {
                        overlap = true;
                    }
                    else
                    {
                        forAll(bubbleCenters, i)
                        {
                            scalar distance = mag(bubbleCenters[i]-newBubble);

                            if (distance < 2.0*bubbleRadius + minBubbleSpacing)
                            {
                                overlap = true;
                                break;
                            }
                        }
                    }

                    // Don't initialise bubbles right below cylinder,
                    // prevents break-up at startup.
                    if
                    (
                           newBubble.z() < 0.05
                        && newBubble.z() > 0.02
                        && Foam::mag(newBubble.x()) < 0.006
                    )
                    {
                        overlap = true;
                    }

                    if (!overlap)
                    {
                        bubbleCenters.append(newBubble);
                    }
                }

                iter++;

                if (iter >= 10000*nBubbles)
                {
                    FatalError
                        << "Maximum iterations exceeded: " << iter
                        << endl;
                    FatalError.exit();
                }
            }
        }

        Pstream::template scatter<vectorList>(bubbleCenters);

        forAllCells(alpha, i, j, k)
            for (int vi = 0; vi < 8; vi++)
                forAll(bubbleCenters, b)
                    alpha(i,j,k) +=
                        0.125
                        *
                        (mag(v(i,j,k)[vi] - bubbleCenters[b]) <= bubbleRadius);

        alpha.correctBoundaryConditions();
    }

    return true;
}

bool PILULE::execute()
{
    const bool staggeredSolver =
        runTime_.foundObject<staggeredScalarField>("U");

    const colocatedVectorDirection& Uc =
        runTime_.lookupObject<colocatedVectorField>
        (
            staggeredSolver ? "Uc" : "U"
        )[0][0];

    colocatedScalarDirection& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

    const colocatedScalarDirection& cv =
        fvMsh_.metrics<colocated>().cellVolumes()[0][0];

    const colocatedScalarDirection& ibmMask =
        fvMsh_.metrics<colocated>().mask()[0][0];

    // Clear alpha stuck in IBM

    const objectRegistry& db = fvMsh_.db();

    label index = 0;
    word name("alpha.0");

    while (db.foundObject<colocatedScalarField>(name))
    {
        colocatedScalarField& alphai =
            db.lookupObjectRef<colocatedScalarField>(name);

        forAllCells(alphai,i,j,k)
        {
            if (ibmMask(i,j,k))
            {
                removedVol_ += alphai(i,j,k) * cv(i,j,k);

                volLossIBM_(i,j,k) += alphai(i,j,k);

                alphai(i,j,k) = Zero;
            }
        }

        name = IOobject::groupName("alpha", Foam::name(++index));
    }

    forAllCells(alpha,i,j,k)
    {
        if (ibmMask(i,j,k))
        {
            alpha(i,j,k) = Zero;
        }
    }

    U1_ = Uc * alpha;
    Usqr1_ = cmptSqr(Uc) * alpha;

    U2_ = Uc * (1.0 - alpha);
    Usqr2_ = cmptSqr(Uc) * (1.0 - alpha);

    const scalar UyInt1
    (
        gSum(Uc.component(vector::Z)*cv*alpha*(1.0-ibmMask))
    );

    const scalar UyInt2
    (
        gSum(Uc.component(vector::Z)*cv*(1.0-alpha)*(1.0-ibmMask))
    );

    const scalar J1 = UyInt1/fluidVol_;
    const scalar J2 = UyInt2/fluidVol_;

    const scalar alphaH = J1/Foam::max(J1+J2,1e-10);

    // PI controller for source term to get desired J2 (J_liquid)
    if (runTime_.time().value() > 1.0) // PI controller startTime
    {
        scalar error = targetJ2_ - J2;

        scalar adjustment =
            kP_ * (error - oldError_) + kI_ * error
          * runTime_.time().deltaTValue();

        forcing_ += adjustment;
        reduce(forcing_, minOp<scalar>());

        oldError_ = error;

        const vector fv(0,0,forcing_);

        // Update the global forcing in fluid cells

        if (staggeredSolver)
        {
            const tensor base =
                fvMsh_.msh().cast<rectilinearMesh>().base();

            staggeredScalarField& exSource =
                runTime_.lookupObjectRef<staggeredScalarField>("source");

            exSource =
                fvMsh_.metrics<staggered>().fluidMask() * list(base & fv);
        }
        else
        {
            colocatedVectorField& exSource =
                runTime_.lookupObjectRef<colocatedVectorField>("source");

            exSource = fvMsh_.metrics<colocated>().fluidMask() * fv;
        }
    }

    if (Pstream::master())
    {
        filePtr_()
                << runTime_.time().value()
                << "; " << J1 << "; " << J2
                << "; " << alphaH  << "; " << forcing_ << nl;
    }

    piluleDict_.set("alphaH", alphaH);

    Info
        << "J_gas: " << J1
        << ", J_liq: " << J2
        << ", alpha_H: " << alphaH
        << ", Removed volume: " << returnReduce(removedVol_, sumOp<scalar>())
        << endl;

    return true;
}

}

}

}

}
