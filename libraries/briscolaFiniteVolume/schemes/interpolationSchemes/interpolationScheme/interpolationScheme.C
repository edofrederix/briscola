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
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("interpolationSchemes").subDict(name)
    );

    const word interpolationSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(interpolationSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation scheme " << interpolationSchemeType
            << nl << nl << "Valid interpolation schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<interpolationScheme<Type,MeshType>>
    (
        cstrIter()(dict, fvMsh)
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
            0.5*(field(ijk) + field(ijkm));
    }

    return tInterp;
}

template<class Type>
tmp<meshField<Type,staggered>> stagInterp
(
    const meshField<FaceSpace<Type>,colocated>& field
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

    forAllCells(field, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector ijk0(ijk + units[0]);
        const labelVector ijk1(ijk + units[1]);
        const labelVector ijk2(ijk + units[2]);

        Interp(0,ijk) = field(ijk).left();
        Interp(0,ijk0) = field(ijk).right();

        Interp(1,ijk) = field(ijk).bottom();
        Interp(1,ijk1) = field(ijk).top();

        Interp(2,ijk) = field(ijk).aft();
        Interp(2,ijk2) = field(ijk).fore();
    }

    return tInterp;
}

template<class Type>
tmp<meshField<Type,staggered>> stagVectorInterp
(
    const meshField<Vector<Type>,colocated>& field
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
            0.5*(field(ijk)[d] + field(ijkm)[d]);
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
