#include "linearRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void linearRestrictionScheme<Type,MeshType>::setWeights()
{
    for (label l = 0; l < this->fvMsh().size()-1; l++)
    {
        const labelVector R(this->fvMsh()[l+1].R());

        for (label d = 0; d < MeshType::numberOfDirections; d++)
        {
            const meshDirection<vector,MeshType>& ccf =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l][d];

            const meshDirection<vector,MeshType>& ccc =
                this->fvMsh().template metrics<MeshType>()
               .cellCenters()[l+1][d];

            meshDirection<vertexScalar,MeshType>& weights = weights_[l+1][d];

            weights = Zero;

            forAllCells(ccc, i, j, k)
            {
                vertexVector vertices;

                int q = 0;
                for (int c = 0; c < 2; c++)
                    for (int b = 0; b < 2; b++)
                        for (int a = 0; a < 2; a++)
                            vertices[q++] =
                                ccf
                                (
                                    i*R.x() + a,
                                    j*R.y() + b,
                                    k*R.z() + c
                                );

                // Calculate the tri-linear interpolation point

                const vector w = interpolationWeights(ccc(i,j,k), vertices);

                // Set the weights

                weights(i,j,k) =
                    vertexScalar
                    (
                        (1.0-w.x()) * (1.0-w.y()) * (1.0-w.z()),
                        (    w.x()) * (1.0-w.y()) * (1.0-w.z()),
                        (1.0-w.x()) * (    w.y()) * (1.0-w.z()),
                        (    w.x()) * (    w.y()) * (1.0-w.z()),
                        (1.0-w.x()) * (1.0-w.y()) * (    w.z()),
                        (    w.x()) * (1.0-w.y()) * (    w.z()),
                        (1.0-w.x()) * (    w.y()) * (    w.z()),
                        (    w.x()) * (    w.y()) * (    w.z())
                    );

                // Important to remove small weights to avoid using
                // uninitialized memory!

                for (int vi = 0; vi < 8; vi++)
                    if (weights(i,j,k)[vi] <= 1e-12)
                        weights(i,j,k)[vi] = 0.0;
            }
        }
    }
}

template<class Type, class MeshType>
linearRestrictionScheme<Type,MeshType>::linearRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<Type,MeshType>(fvMsh, is),
    weights_
    (
        "weights",
        this->fvMsh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    )
{
    setWeights();
}

template<class Type, class MeshType>
void linearRestrictionScheme<Type,MeshType>::restrict
(
    meshLevel<Type,MeshType>& coarse,
    const meshLevel<Type,MeshType>& fine,
    const bool scale
)
{
    const labelVector R(coarse.lvl().R());

    meshLevel<vertexScalar,MeshType>& weights =
        weights_[coarse.levelNum()];

    const_cast<meshLevel<Type,MeshType>&>(fine).correctAggData();

    if (scale)
    {
        const meshLevel<scalar,MeshType>& cvc =
            this->fvMsh().template
            metrics<MeshType>().cellVolumes()[coarse.levelNum()];

        const meshLevel<scalar,MeshType>& icvf =
            this->fvMsh().template
            metrics<MeshType>().inverseCellVolumes()[fine.levelNum()];

        forAllCells(coarse, d, i, j, k)
        {
            coarse(d,i,j,k) = Zero;

            const label ox = R.x() > 1;
            const label oy = R.y() > 1;
            const label oz = R.z() > 1;

            int q = 0;
            for (int c = 0; c < 2; c++)
                for (int b = 0; b < 2; b++)
                    for (int a = 0; a < 2; a++)
                        if (weights(d,i,j,k)[q++] != 0.0)
                            coarse(d,i,j,k) +=
                                weights(d,i,j,k)[q-1]
                              * fine
                                (
                                    d,
                                    i*R.x() + a*ox,
                                    j*R.y() + b*oy,
                                    k*R.z() + c*oz
                                )
                              * icvf
                                (
                                    d,
                                    i*R.x() + a*ox,
                                    j*R.y() + b*oy,
                                    k*R.z() + c*oz
                                );

            coarse(d,i,j,k) *= cvc(d,i,j,k);
        }
    }
    else
    {
        forAllCells(coarse, d, i, j, k)
        {
            coarse(d,i,j,k) = Zero;

            const label ox = R.x() > 1;
            const label oy = R.y() > 1;
            const label oz = R.z() > 1;

            int q = 0;
            for (int c = 0; c < 2; c++)
                for (int b = 0; b < 2; b++)
                    for (int a = 0; a < 2; a++)
                        if (weights(d,i,j,k)[q++] != 0.0)
                            coarse(d,i,j,k) +=
                                weights(d,i,j,k)[q-1]
                              * fine
                                (
                                    d,
                                    i*R.x() + a*ox,
                                    j*R.y() + b*oy,
                                    k*R.z() + c*oz
                                );
        }
    }
}

}

}

}
