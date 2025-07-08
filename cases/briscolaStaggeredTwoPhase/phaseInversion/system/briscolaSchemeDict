FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaSchemesDict;
}

laplacianSchemes
{
    default linearGauss;
}

divergenceSchemes
{
    div(phi,U) limitedGauss twoPhase vanLeer;
    // default midPoint
}

gradientSchemes
{
    default midPointGauss;
}

faceGradientSchemes
{
    default linear;
}

faceFluxSchemes
{
    default midPoint;
}

interpolationSchemes
{
    default midPoint;
}

reconstructionSchemes
{
    default midPoint;
}

rkScheme
{
    scheme Ascher222;
}