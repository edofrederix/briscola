FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      briscolaSchemesDict;
}

ddtSchemes
{
    default Euler;
}

laplacianSchemes
{
    default VARLSCHEME;
}

divergenceSchemes
{
    default VARDSCHEME;
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

rkScheme
{
    scheme VARRKSCHEME;
}
