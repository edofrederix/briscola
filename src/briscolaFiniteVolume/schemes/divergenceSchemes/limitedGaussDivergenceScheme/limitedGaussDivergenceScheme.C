#include "limitedGaussDivergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
limitedGaussDivergenceScheme<Type,MeshType>::limitedGaussDivergenceScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    divergenceScheme<stencil,Type,MeshType>(fvMsh, is),
    limiter_(limiterScheme<Type,MeshType>::New(fvMsh,is))
{}

template<class Type, class MeshType>
tmp<linearSystem<stencil,Type,MeshType>>
limitedGaussDivergenceScheme<Type,MeshType>::imDiv
(
    const meshField<faceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    phi.restrict();
    field.restrict();

    tmp<linearSystem<stencil,Type,MeshType>> tSys
    (
        new linearSystem<stencil,Type,MeshType>
        (
            word("div(" + phi.name() + "," + field.name() + ")"),
            const_cast<meshField<Type,MeshType>&>(field)
        )
    );

    linearSystem<stencil,Type,MeshType>& sys = tSys.ref();

    meshField<stencil,MeshType>& A = sys.A();

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    const meshField<faceScalar,MeshType> psi(limiter_->psi(phi,field));

    forAllCells(A, l, d, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        A(l,d,i,j,k)[0] = Zero;
        #endif

        for (label f = 0; f < 6; f++)
        {
            A(l,d,i,j,k)[f+1] =
                factor*phi(l,d,i,j,k)[f]
              * (
                    psi(l,d,i,j,k)[f]*fwn(l,d,i,j,k)[f]
                  + (1.0 - psi(l,d,i,j,k)[f])*(phi(l,d,i,j,k)[f] < 0)
                );

            A(l,d,i,j,k)[0] +=
                factor*phi(l,d,i,j,k)[f]
              * (
                    psi(l,d,i,j,k)[f]*fwc(l,d,i,j,k)[f]
                  + (1.0 - psi(l,d,i,j,k)[f])*(phi(l,d,i,j,k)[f] >= 0)
                );
        }
    }

    #ifdef NO_BLOCK_ZERO_INIT

    meshField<Type,MeshType>& b = sys.b();
    forAllCells(b, l, d, i, j, k)
        b(l,d,i,j,k) = Zero;

    #endif

    phi.makeShallow();
    field.makeShallow();

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
limitedGaussDivergenceScheme<Type,MeshType>::exDiv
(
    const meshField<faceScalar,MeshType>& phi,
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tDiv
    (
        new meshField<Type,MeshType>
        (
            "div("+phi.name()+","+field.name()+")",
            phi.fvMsh()
        )
    );

    meshField<Type,MeshType>& Div = tDiv.ref();

    const meshField<faceScalar,MeshType>& fwc =
        field.fvMsh().template metrics<MeshType>().faceWeightsCenter();

    const meshField<faceScalar,MeshType>& fwn =
        field.fvMsh().template metrics<MeshType>().faceWeightsNeighbor();

    const meshField<scalar,MeshType>& cv =
        phi.fvMsh().template metrics<MeshType>().cellVolumes();

    const meshField<faceScalar,MeshType> psi(limiter_->psi(phi,field));

    forAllCells(Div, d, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Div(d,i,j,k) = Zero;
        #endif

        for (label f = 0; f < 6; f++)
            Div(d,i,j,k) +=
                phi(d,i,j,k)[f]
              * (
                    field(d,i,j,k)
                  * (
                        psi(d,i,j,k)[f]*fwc(d,i,j,k)[f]
                      + (1.0 - psi(d,i,j,k)[f])*(phi(d,i,j,k)[f] >= 0)
                    )
                  + field(d,neighbor(i,j,k,f))
                  * (
                        psi(d,i,j,k)[f]*fwn(d,i,j,k)[f]
                      + (1.0 - psi(d,i,j,k)[f])*(phi(d,i,j,k)[f] < 0)
                    )
                )
              / cv(d,i,j,k);
    }

    return tDiv;
}

}

}

}
