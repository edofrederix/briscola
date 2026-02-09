#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"
#include "IO.H"

#include "restrictionSchemes.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class MeshType>
void testAggMetrics(const fvMesh& fvMsh)
{
    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    const faceField<scalar,MeshType>& delta =
        fvMsh.metrics<MeshType>().faceDeltas();

    const tensor base(fvMsh.msh().cast<rectilinearMesh>().base());

    forAll(cc, l)
    {
        const decomposition& decomp = fvMsh[l].decomp();

        if (decomp.aggParent() && fvMsh[l].child().decomp().aggMaster())
        {
            forAll(cc[l], d)
            {
                const vectorBlock& B = cc[l][d].B();

                const scalarBlock& Dx = delta[0][l][d].B();
                const scalarBlock& Dy = delta[1][l][d].B();
                const scalarBlock& Dz = delta[2][l][d].B();

                for (int i = 1; i < B.l()-1; i++)
                for (int j = 1; j < B.m()-1; j++)
                for (int k = 1; k < B.n()-1; k++)
                {
                    // Cell centers must be increasing along the base directions

                    if (((B(i,j,k) - B(i-1,j,k)) & base.x()) <= 0)
                        FatalErrorInFunction
                            << "Test 1a failed" << abort(FatalError);

                    if (((B(i,j,k) - B(i,j-1,k)) & base.y()) <= 0)
                        FatalErrorInFunction
                            << "Test 1b failed" << abort(FatalError);

                    if (((B(i,j,k) - B(i,j,k-1)) & base.z()) <= 0)
                        FatalErrorInFunction
                            << "Test 1c failed" << abort(FatalError);

                    // Cell centered distances must be equal

                    if (i > 1)
                        if
                        (
                            Foam::mag
                            (
                                Foam::mag(B(i,j,k) - B(i-1,j,k))
                              - Foam::mag(B(i-1,j,k) - B(i-2,j,k))
                            ) > 1e-8
                        )
                            FatalErrorInFunction
                                << "Test 2a failed" << abort(FatalError);

                    if (j > 1)
                        if
                        (
                            Foam::mag
                            (
                                Foam::mag(B(i,j,k) - B(i,j-1,k))
                              - Foam::mag(B(i,j-1,k) - B(i,j-2,k))
                            ) > 1e-8
                        )
                            FatalErrorInFunction
                                << "Test 2b failed" << abort(FatalError);

                    if (k > 1)
                        if
                        (
                            Foam::mag
                            (
                                Foam::mag(B(i,j,k) - B(i,j,k-1))
                              - Foam::mag(B(i,j,k-1) - B(i,j,k-2))
                            ) > 1e-8
                        )
                            FatalErrorInFunction
                                << "Test 2c failed" << abort(FatalError);

                    // Check delta values

                    if
                    (
                        Foam::mag
                        (
                            ((B(i,j,k) - B(i-1,j,k)) & base.x())
                          - 1.0/Dx(i,j,k)
                        ) > 1e-8
                    )
                        FatalErrorInFunction
                            << "Test 3a failed" << abort(FatalError);

                    if
                    (
                        Foam::mag
                        (
                            ((B(i,j,k) - B(i,j-1,k)) & base.y())
                          - 1.0/Dy(i,j,k)
                        ) > 1e-8
                    )
                        FatalErrorInFunction
                            << "Test 3b failed" << abort(FatalError);

                    if
                    (
                        Foam::mag
                        (
                            ((B(i,j,k) - B(i,j,k-1)) & base.z())
                          - 1.0/Dz(i,j,k)
                        ) > 1e-8
                    )
                        FatalErrorInFunction
                            << "Test 3c failed" << abort(FatalError);
                }
            }
        }
    }
}

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

    #include "createBriscolaIO.H"

    // Test cell centers

    testAggMetrics<colocated>(fvMsh);
    testAggMetrics<staggered>(fvMsh);
}
