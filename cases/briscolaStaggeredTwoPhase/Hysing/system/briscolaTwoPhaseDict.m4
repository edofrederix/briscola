FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaTwoPhaseDict;
}

type    twoPhaseVof;

g       (0 -0.98 0);

rho1    1000.0;
rho2    VARRHO2;
mu1     10.0;
mu2     VARMU2;

normalScheme
{
    type    MYC;
}

surfaceTensionScheme
{
    type    Brackbill;

    curvatureScheme
    {
        type    SHF;
    }

    sigma   VARSIGMA;
}

vof
{
    type    splitAdvection;
}

reduced true;
