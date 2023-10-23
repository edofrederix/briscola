#include "LVE.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

const scalar LVE::tol_ = 1e-8;
const label LVE::maxIter_ = 100;

scalar LVE::fluxVolumeLVE
(
    const vertexVector& v,
    const scalar& Vtotal,
    const scalar& Vf,
    const label& u,
    const label& l
) const
{
    vectorList x(8);
    vectorList e(8);

    for (int iv = 0; iv < 4; iv++)
    {
        x[vertexNumsInFace[u][iv]] = v[vertexNumsInFace[u][iv]];
        x[vertexNumsInFace[l][iv]] = v[vertexNumsInFace[l][iv]];

        e[vertexNumsInFace[u][iv]] = v[vertexNumsInFace[l][iv]] - v[vertexNumsInFace[u][iv]];
        e[vertexNumsInFace[l][iv]] = Zero;
    }

    scalar alpha0 = 0;
    scalar alpha1 = 0;
    scalar alpha2 = 0;
    scalar alpha3 = 0;

    for (int i = 0; i < numberOfTets; i++)
    {
        vectorList x1(3);
        vectorList e1(3);

        for (int j = 0; j < 3; j++)
        {
            x1[j] = x[tetDecomp[i][j]] - x[tetDecomp[i][3]];
            e1[j] = e[tetDecomp[i][j]] - e[tetDecomp[i][3]];
        }

        alpha0 += (x1[0] & (x1[1] ^ x1[2]));
        alpha1 += (e1[0] & (x1[1] ^ x1[2]))
                + (x1[0] & (e1[1] ^ x1[2]))
                + (x1[0] & (x1[1] ^ e1[2]));
        alpha2 += (x1[0] & (e1[1] ^ e1[2]))
                + (e1[0] & (x1[1] ^ e1[2]))
                + (e1[0] & (e1[1] ^ x1[2]));
        alpha3 += (e1[0] & (e1[1] ^ e1[2]));

    }

    scalar C = - LVE::exactCubicSolver
    (
        Vf,
        alpha0,
        alpha1,
        alpha2,
        alpha3,
        0,
        1,
        0,
        Vtotal,
        Vtotal
    );

    if (C > 1)
    {
        C = 1;
    }
    else if (C < 0)
    {
        C = 0;
    }

    return C;

}

}

}

}


