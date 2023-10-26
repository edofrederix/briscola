#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
interpolationScheme<Type,MeshType>::interpolationScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class Type, class MeshType>
interpolationScheme<Type,MeshType>::interpolationScheme
(
    const interpolationScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
interpolationScheme<Type,MeshType>::~interpolationScheme()
{}

template<class Type, class MeshType>
autoPtr<interpolationScheme<Type,MeshType>>
interpolationScheme<Type,MeshType>::New
(
    const word name,
    const fvMesh& fvMsh
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation scheme " << name << nl << nl
            << "Valid interpolation schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<interpolationScheme<Type,MeshType>>
    (
        cstrIter()(dictionary(), fvMsh)
    );
}

template<class Type>
tmp<meshField<Type,staggered>> stagInterp
(
    const meshField<Type,colocated>& field
)
{
    tmp<meshField<Type,staggered>> tInterp
    (
        new meshField<Type,staggered>
        (
            "stagInterp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,staggered>& Interp = tInterp.ref();

    Interp = Zero;

    forAllDirections(Interp, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector ijkm(ijk - units[d]);

        // Staggered cell center is exactly in the middle of the two colocated
        // cell centers

        Interp(d,ijk) =
            0.5*(field(d,ijk) + field(d,ijkm));
    }

    return tInterp;
}

template<class Type>
tmp<meshField<FaceSpace<Type>,staggered>> stagFaceInterp
(
    const meshField<Type,colocated>& field
)
{
    tmp<meshField<FaceSpace<Type>,staggered>> tInterp
    (
        new meshField<FaceSpace<Type>,staggered>
        (
            "stagFaceInterp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<FaceSpace<Type>,staggered>& Interp = tInterp.ref();

    Interp = Zero;

    forAllDirections(Interp, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector ijkm(ijk - units[d]);

        for (int f = 0; f < 6; f++)
        {
            if (f == d*2)
            {
                // Staggered face coincides with lower colocated center

                Interp(d,ijk)[f] = field(0,ijkm);
            }
            else if (f == d*2+1)
            {
                // Staggered face coincides with upper colocated center

                Interp(d,ijk)[f] = field(0,ijk);
            }
            else if (f%2 == 0)
            {
                // Even face number

                Interp(d,ijk)[f] =
                    0.25
                  * (
                        field(0,ijk)
                      + field(0,ijkm)
                      + field(0,ijk  - units[f/2])
                      + field(0,ijkm - units[f/2])
                    );
            }
            else
            {
                // Odd face number

                Interp(d,ijk)[f] =
                    0.25
                  * (
                        field(0,ijk)
                      + field(0,ijkm)
                      + field(0,ijk  + units[f/2])
                      + field(0,ijkm + units[f/2])
                    );
            }
        }
    }

    return tInterp;
}

}

}

}
