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
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
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
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "interpolationSchemes", schemeName);

    word interpolationSchemeType;
    is >> interpolationSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(interpolationSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation scheme " << interpolationSchemeType
            << nl << nl << "Valid interpolation schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<interpolationScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, is)
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

    forAllCells(Interp, d, i, j, k)
        Interp(d,i,j,k) =
            0.5*(field(i,j,k) + field(lowerNeighbor(i,j,k,d)));

    return tInterp;
}

template<class Type>
tmp<meshField<Type,staggered>> stagInterp
(
    const tmp<meshField<Type,colocated>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<Type,staggered>> tInterp
    (
        stagInterp(tField())
    );

    if (tField.isTmp())
        tField.clear();

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

    forAllFaces(Interp, d, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,d));

        if (d == fd)
        {
            Interp(d,ijk)[fd*2] = field(nei);
        }
        else
        {
            Interp(d,ijk)[fd*2] =
                0.25
              * (
                    field(ijk)
                  + field(lowerNeighbor(ijk,fd))
                  + field(nei)
                  + field(lowerNeighbor(nei,fd))
                );
        }

        Interp(d,lowerNeighbor(i,j,k,fd))[fd*2+1] = Interp(d,ijk)[fd*2];
    }

    return tInterp;
}

template<class Type>
tmp<meshField<FaceSpace<Type>,staggered>> stagFaceInterp
(
    const tmp<meshField<Type,colocated>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<FaceSpace<Type>,staggered>> tInterp
    (
        stagFaceInterp(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}

template<class Type>
tmp<meshField<Type,colocated>> coloInterp
(
    const meshField<Type,staggered>& field
)
{
    tmp<meshField<Type,colocated>> tInterp
    (
        new meshField<Type,colocated>
        (
            "coloInterp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,colocated>& Interp = tInterp.ref();

    forAllCells(Interp, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Interp(i,j,k) = Zero;
        #endif

        for (int d = 0; d < 3; d++)
            Interp(i,j,k) +=
                0.5/3.0*(field(d,i,j,k) + field(upperNeighbor(i,j,k,d)));
    }

    return tInterp;
}

template<class Type>
tmp<meshField<Type,colocated>> coloInterp
(
    const tmp<meshField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<Type,colocated>> tInterp
    (
        stagInterp(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}


template<class Type>
tmp<meshField<FaceSpace<Type>,colocated>> coloFaceInterp
(
    const meshField<Type,staggered>& field
)
{
    tmp<meshField<FaceSpace<Type>,colocated>> tInterp
    (
        new meshField<FaceSpace<Type>,colocated>
        (
            "coloFaceInterp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<FaceSpace<Type>,colocated>& Interp = tInterp.ref();

    forAllFaces(Interp, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        Interp(ijk)[fd*2  ] = field(fd,i,j,k);
        Interp(nei)[fd*2+1] = Interp(ijk)[fd*2];
    }

    return tInterp;
}

template<class Type>
tmp<meshField<FaceSpace<Type>,colocated>> coloFaceInterp
(
    const tmp<meshField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<FaceSpace<Type>,colocated>> tInterp
    (
        coloFaceInterp(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tInterp;
}

}

}

}
