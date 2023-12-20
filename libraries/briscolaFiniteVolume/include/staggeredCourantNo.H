scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const staggeredScalarLevel& cvl =
        fvMsh.metrics<staggered>().cellVolumes().mshPart();

    staggeredScalarLevel sumPhi(fvMsh,0,Zero);

    forAllFaces(sumPhi, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        sumPhi(d,ijk) += mag(phi(d,ijk)[fd]);
        sumPhi(d,nei) += mag(phi(d,ijk)[fd]);
    }

    CoNum = max(0.5*gMax(sumPhi/cvl)*runTime.deltaTValue());

    meanCoNum =
        max(0.5*(gSum(sumPhi)/gSum(cvl))*runTime.deltaTValue());
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;