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
    default VARDSCHEME;
}

gradientSchemes
{
    default VARGRADSCHEME;
}

faceGradientSchemes
{
    default linear;
}

faceFluxSchemes
{
    default linear;
}

rkScheme
{
    scheme VARRKSCHEME;
}
