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

        forAll(sizes, proc)
        {
            const label size = sizes[proc];

            buffer.clear();
            buffer.setSize(size);

            forAll(buffer, i)
            {
                file>> buffer[i];
            }

            if (proc == Pstream::myProcNo())
            {
                data = buffer;
            }
        }

        nextLine(file);
    }
    else
    {
        label offset = 0;

        for(int i = 0; i < Pstream::myProcNo(); i++)
            offset += sizes[i];

        file.ignore(offset*sizeof(floatScalar));

        file.read
        (
            reinterpret_cast<char*>(data.begin()),
            data.byteSize()
        );

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
}

template<class Type, class MeshType>
void IO::readScalarField
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    meshDirection<Type,MeshType>& D
) const
{
    List<floatScalar> data(D.size());

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

    forAllCells(D, i, j, k)
        D(i,j,k) = data[c++];
}

template<class Type, class MeshType>
void IO::readArrayField
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nComponents);

    List<floatScalar> data(D.size()*n);

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

    forAllCells(D, i, j, k)
        for (int ii = 0; ii < n; ii++)
            D(i,j,k)[ii] = data[c++];
}

template<class Type, class MeshType>
void IO::readArrayArrayField
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nComponents);
    const label m(pTraits<typename Type::cmpt>::nComponents);

    List<floatScalar> data(D.size()*m*n);

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

    forAllCells(D, i, j, k)
        for (int ii = 0; ii < n; ii++)
            for (int jj = 0; jj < m; jj++)
                D(i,j,k)[ii][jj] = data[c++];
}

// Instantiate

#define READTYPEFIELD(FUNC,TYPE,MESHTYPE)                                       \
                                                                                \
template<>                                                                      \
void IO::readField                                                              \
(                                                                               \
    autoPtr<std::ifstream>& filePtr,                                            \
    const bool ascii,                                                           \
    meshDirection<TYPE,MESHTYPE>& D                                             \
) const                                                                         \
{                                                                               \
    FUNC(filePtr,ascii,D);                                                      \
}                                                                               \
                                                                                \
template void IO::FUNC                                                          \
(                                                                               \
    autoPtr<std::ifstream>&,                                                    \
    const bool ascii,                                                           \
    meshDirection<TYPE,MESHTYPE>&                                               \
) const;

READTYPEFIELD(readScalarField,scalar,colocated)
READTYPEFIELD(readScalarField,label,colocated)
READTYPEFIELD(readArrayField,vector,colocated)
READTYPEFIELD(readArrayField,tensor,colocated)
READTYPEFIELD(readArrayField,diagTensor,colocated)
READTYPEFIELD(readArrayField,symmTensor,colocated)
READTYPEFIELD(readArrayField,sphericalTensor,colocated)
READTYPEFIELD(readArrayField,faceScalar,colocated)
READTYPEFIELD(readArrayField,edgeScalar,colocated)
READTYPEFIELD(readArrayField,vertexScalar,colocated)
READTYPEFIELD(readArrayArrayField,faceVector,colocated)
READTYPEFIELD(readArrayArrayField,edgeVector,colocated)
READTYPEFIELD(readArrayArrayField,vertexVector,colocated)

READTYPEFIELD(readScalarField,scalar,staggered)
READTYPEFIELD(readScalarField,label,staggered)
READTYPEFIELD(readArrayField,vector,staggered)
READTYPEFIELD(readArrayField,tensor,staggered)
READTYPEFIELD(readArrayField,diagTensor,staggered)
READTYPEFIELD(readArrayField,symmTensor,staggered)
READTYPEFIELD(readArrayField,sphericalTensor,staggered)
READTYPEFIELD(readArrayField,faceScalar,staggered)
READTYPEFIELD(readArrayField,edgeScalar,staggered)
READTYPEFIELD(readArrayField,vertexScalar,staggered)
READTYPEFIELD(readArrayArrayField,faceVector,staggered)
READTYPEFIELD(readArrayArrayField,edgeVector,staggered)
READTYPEFIELD(readArrayArrayField,vertexVector,staggered)

#undef READTYPEFIELD

}

}

}