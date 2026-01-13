#include "IO.H"
#include "colocatedFields.H"
#include "staggeredFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Not a templated function but used mostly by templated functions

void IO::readList
(
    autoPtr<std::ifstream>& filePtr,
    List<floatScalar>& data,
    const bool ascii,
    const label tag
) const
{
    labelList sizes(Pstream::nProcs());
    sizes[Pstream::myProcNo()] = data.size();
    Pstream::gatherList(sizes);
    Pstream::scatterList(sizes);

    std::ifstream& file = filePtr();

    if (ascii)
    {
        List<floatScalar> buffer;

        forAll(sizes, p)
        if (!partitioned_ || p == Pstream::myProcNo())
        {
            const label size = sizes[p];

            buffer.clear();
            buffer.setSize(size);

            forAll(buffer, i)
            {
                file>> buffer[i];
            }

            if (p == Pstream::myProcNo())
            {
                data = buffer;
            }
        }

        nextLine(file);
    }
    else
    {
        label offset = 0;

        for(int i = 0; i < (partitioned_ ? 0 : Pstream::myProcNo()); i++)
            offset += sizes[i];

        file.ignore(offset*sizeof(floatScalar));

        file.read
        (
            reinterpret_cast<char*>(data.begin()),
            data.byteSize()
        );

        if (!partitioned_)
            file.ignore((sum(sizes)-offset-data.size())*sizeof(floatScalar));

        nextLine(file);

        #ifdef LITTLEENDIAN
        swapWords
        (
            data.byteSize()/sizeof(label),
            reinterpret_cast<label*>(data.begin())
        );
        #endif
    }

    if (Pstream::parRun())
        returnReduce(1.0,sumOp<scalar>());
}

template<class Type, class MeshType>
void IO::readScalarField
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    meshDirection<Type,MeshType>& D
) const
{
    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_)*D.ghosts;
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_)*D.ghosts;
    const labelVector N = E - S;

    List<floatScalar> data(nStructured(N));

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::readList
    (
        filePtr,
        data,
        ascii,
        tag
    );

    label c = 0;

    for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (structured(i-S.x(),j-S.y(),k-S.z(),N))
                    D(i,j,k) = data[c++];
}

template<class Type, class MeshType>
void IO::readVectorSpaceField
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nComponents);

    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_)*D.ghosts;
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_)*D.ghosts;
    const labelVector N = E - S;

    List<floatScalar> data(nStructured(N)*n);

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::readList
    (
        filePtr,
        data,
        ascii,
        tag
    );

    label c = 0;

    for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (structured(i-S.x(),j-S.y(),k-S.z(),N))
                    for (int ii = 0; ii < n; ii++)
                        setComponent(D(i,j,k),ii) = data[c++];
}

template<class Type, class MeshType>
void IO::readCellSpaceField
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nCsComponents);
    const label m(Type::nComponents);

    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_)*D.ghosts;
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_)*D.ghosts;
    const labelVector N = E - S;

    List<floatScalar> data(nStructured(N)*m*n);

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::readList
    (
        filePtr,
        data,
        ascii,
        tag
    );

    label c = 0;

    for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (structured(i-S.x(),j-S.y(),k-S.z(),N))
                    for (int ii = 0; ii < n; ii++)
                        for (int jj = 0; jj < m; jj++)
                            setComponent(D(i,j,k)[ii],jj) = data[c++];
}

// Instantiate

#define READTYPEFIELD(FUNC,TYPE,MESHTYPE)                                      \
                                                                               \
template<>                                                                     \
void IO::readField                                                             \
(                                                                              \
    autoPtr<std::ifstream>& filePtr,                                           \
    const bool ascii,                                                          \
    meshDirection<TYPE,MESHTYPE>& D                                            \
) const                                                                        \
{                                                                              \
    FUNC(filePtr,ascii,D);                                                     \
}                                                                              \
                                                                               \
template void IO::FUNC                                                         \
(                                                                              \
    autoPtr<std::ifstream>&,                                                   \
    const bool ascii,                                                          \
    meshDirection<TYPE,MESHTYPE>&                                              \
) const;

READTYPEFIELD(readScalarField,scalar,colocated)
READTYPEFIELD(readScalarField,label,colocated)
READTYPEFIELD(readVectorSpaceField,vector,colocated)
READTYPEFIELD(readVectorSpaceField,tensor,colocated)
READTYPEFIELD(readVectorSpaceField,diagTensor,colocated)
READTYPEFIELD(readVectorSpaceField,symmTensor,colocated)
READTYPEFIELD(readVectorSpaceField,sphericalTensor,colocated)
READTYPEFIELD(readCellSpaceField,faceScalar,colocated)
READTYPEFIELD(readCellSpaceField,edgeScalar,colocated)
READTYPEFIELD(readCellSpaceField,vertexScalar,colocated)
READTYPEFIELD(readCellSpaceField,faceVector,colocated)
READTYPEFIELD(readCellSpaceField,edgeVector,colocated)
READTYPEFIELD(readCellSpaceField,vertexVector,colocated)

READTYPEFIELD(readScalarField,scalar,staggered)
READTYPEFIELD(readScalarField,label,staggered)
READTYPEFIELD(readVectorSpaceField,vector,staggered)
READTYPEFIELD(readVectorSpaceField,tensor,staggered)
READTYPEFIELD(readVectorSpaceField,diagTensor,staggered)
READTYPEFIELD(readVectorSpaceField,symmTensor,staggered)
READTYPEFIELD(readVectorSpaceField,sphericalTensor,staggered)
READTYPEFIELD(readCellSpaceField,faceScalar,staggered)
READTYPEFIELD(readCellSpaceField,edgeScalar,staggered)
READTYPEFIELD(readCellSpaceField,vertexScalar,staggered)
READTYPEFIELD(readCellSpaceField,faceVector,staggered)
READTYPEFIELD(readCellSpaceField,edgeVector,staggered)
READTYPEFIELD(readCellSpaceField,vertexVector,staggered)

#undef READTYPEFIELD

}

}

}