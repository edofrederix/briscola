scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const staggeredScalarLevel& cvl =
        fvMsh.metrics<staggered>().cellVolumes().level();

    staggeredScalarLevel sumPhi(fvMsh,0,Zero);

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        sumPhi(d,ijk) += mag(phi[fd](d,ijk));
        sumPhi(d,nei) += mag(phi[fd](d,ijk));
    }

    CoNum = max(0.5*gMax(sumPhi/cvl)*runTime.deltaTValue());

    meanCoNum =
        max(0.5*(gSum(sumPhi)/gSum(cvl))*runTime.deltaTValue());
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;