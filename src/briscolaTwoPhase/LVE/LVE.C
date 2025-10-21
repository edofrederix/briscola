#include "LVE.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

const scalar LVE::tol_ = 1e-12;
const label LVE::maxIter_ = 100;

scalar LVE::fluxVolumeLVE
(
    const vertexVector& v,
    const scalar V,
    const scalar Vf,
    const label u,
    const label l
) const
{
    vectorList x(8);
    vectorList e(8);

    for (int i = 0; i < 4; i++)
    {
        x[vertexNumsInFace[u][i]] = v[vertexNumsInFace[u][i]];
        x[vertexNumsInFace[l][i]] = v[vertexNumsInFace[l][i]];

        e[vertexNumsInFace[u][i]] =
            v[vertexNumsInFace[l][i]]
          - v[vertexNumsInFace[u][i]];

        e[vertexNumsInFace[l][i]] = Zero;
    }

    scalarList alpha(4, 0.0);

    for (int i = 0; i < numberOfTets; i++)
    {
        vectorList w(3);
        vectorList E(3);

        for (int j = 0; j < 3; j++)
        {
            w[j] = x[tetDecomp[i][j]] - x[tetDecomp[i][3]];
            E[j] = E[tetDecomp[i][j]] - E[tetDecomp[i][3]];
        }

        alpha[0] += (w[0] & (w[1] ^ w[2]));
        alpha[3] += (E[0] & (E[1] ^ E[2]));

        alpha[1] +=
            (E[0] & (w[1] ^ w[2]))
          + (w[0] & (E[1] ^ w[2]))
          + (w[0] & (w[1] ^ E[2]));

        alpha[2] +=
            (w[0] & (E[1] ^ E[2]))
          + (E[0] & (w[1] ^ E[2]))
          + (E[0] & (E[1] ^ w[2]));
    }

    const scalar C =
        LVE::exactCubicSolver
        (
            Vf,
            alpha,
            0,
            1,
            0,
            V,
            V
        )
      + 1.0;

    return Foam::max(Foam::min(1.0, C), 0.0);
}

}

}

}
