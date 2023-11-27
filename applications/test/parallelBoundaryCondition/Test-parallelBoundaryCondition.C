#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            runTime.system()/"briscolaMeshDict",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvMesh fvMsh(meshDict, runTime);

    colocatedVectorField f
    (
        "f",
        fvMsh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        true
    );

    const colocatedVectorField& cc =
        fvMsh.metrics<colocated>().cellCenters();

    forAllLevels(f, l, d, i, j, k)
        f(l,d,i,j,k) = cc(l,d,i,j,k);

    f.correctBoundaryConditions();

    forAllLevels(f, l, d, i, j, k)
        if (f(l,d,i,j,k) != cc(l,d,i,j,k))
            FatalError << "test 1 failed" << endl;

    forAll(f.boundaryConditions(), b)
    {
        const boundaryCondition<vector,colocated>& bc =
            f.boundaryConditions()[b];

        if (bc.type() == "parallel")
        {
            const labelVector bo(bc.offset());

            forAll(f, l)
            {
                forAll(f[l], d)
                {
                    const labelVector S(fvMsh.template S<colocated>(l,d,bo));
                    const labelVector E(fvMsh.template E<colocated>(l,d,bo));

                    labelVector ijk;

                    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    {
                        // Domain size is 1x1x1. Uniform parallel partitioning
                        // assumed.

                        const scalar cellSize =
                            1.0/pow(Pstream::nProcs(),1.0/3.0)/f.N(l,d).x();

                        const vector cell = f(l,d,ijk);
                        const vector ghost = f(l,d,ijk+bo);
                        const vector target = cell + vector(bo)*cellSize;

                        if
                        (
                            target.x() > ghost.x()+1e-8
                         || target.x() < ghost.x()-1e-8
                         || target.y() > ghost.y()+1e-8
                         || target.y() < ghost.y()-1e-8
                         || target.z() > ghost.z()+1e-8
                         || target.z() < ghost.z()-1e-8
                        )
                        {
                            FatalError  << "Cell size at this level = "
                                        << cellSize << nl
                                        << "cell = " << cell << nl
                                        << "ghost = " << ghost << nl
                                        << "target = " << target << nl
                                        << "offset = " << bo << nl
                                        << "test 2 failed" << endl;

                            FatalError.exit();
                        }
                    }
                }
            }
        }
    }
}
