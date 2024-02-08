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

    forAllCells(Interp, d, i, j, k)
        Interp(d,i,j,k) =
            0.5*(field(i,j,k) + field(lowerNei(i,j,k,d)));

    return tInterp;
}

template<class Type>
tmp<meshField<LowerFaceSpace<Type>,staggered>> stagFaceInterp
(
    const meshField<Type,colocated>& field
)
{
    tmp<meshField<LowerFaceSpace<Type>,staggered>> tInterp
    (
        new meshField<LowerFaceSpace<Type>,staggered>
        (
            "stagFaceInterp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<LowerFaceSpace<Type>,staggered>& Interp = tInterp.ref();

    Interp = Zero;

    forAllFaces(Interp, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[d]);

        if (d == fd)
        {
            Interp(d,ijk)[fd] = field(nei);
        }
        else
        {
            Interp(d,ijk)[fd] =
                0.25
              * (
                    field(ijk)
                  + field(ijk-units[fd])
                  + field(nei)
                  + field(nei-units[fd])
                );
        }
    }

    return tInterp;
}

}

}

}
