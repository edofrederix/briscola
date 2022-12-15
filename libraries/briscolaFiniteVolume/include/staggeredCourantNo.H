scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const staggeredScalarLevel& cvl =
        fvMsh.metrics<staggered>().cellVolumes()[0];

    staggeredScalarLevel sumPhi(fvMsh,0);

    forAll(sumPhi, d)
    {
        staggeredScalarDirection& sumPhid = sumPhi[d];
        const staggeredFaceScalarDirection& phild = phi[0][d];

        forAllBlock(sumPhid, i, j, k)
        {
            sumPhid(i,j,k) =
                mag(phild(i,j,k).left())
              + mag(phild(i,j,k).right())
              + mag(phild(i,j,k).bottom())
              + mag(phild(i,j,k).top())
              + mag(phild(i,j,k).aft())
              + mag(phild(i,j,k).fore());
        }
    }

    CoNum = max(0.5*gMax(sumPhi/cvl)*runTime.deltaTValue());

    meanCoNum =
        max(0.5*(gSum(sumPhi)/gSum(cvl))*runTime.deltaTValue());
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;