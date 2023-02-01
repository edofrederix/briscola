scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const colocatedFaceScalarDirection& phid = phi[0][0];

    const colocatedScalarDirection& cvd =
        fvMsh.metrics<colocated>().cellVolumes()[0][0];

    colocatedScalarDirection sumPhi(fvMsh,0,0);

    forAllCells(sumPhi, i, j, k)
    {
        sumPhi(i,j,k) =
            mag(phid(i,j,k).left())
          + mag(phid(i,j,k).right())
          + mag(phid(i,j,k).bottom())
          + mag(phid(i,j,k).top())
          + mag(phid(i,j,k).aft())
          + mag(phid(i,j,k).fore());
    }

    CoNum = 0.5*gMax(sumPhi/cvd)*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(cvd))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;