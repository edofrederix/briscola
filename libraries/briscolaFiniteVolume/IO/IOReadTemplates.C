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
    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_);
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_);

    List<floatScalar> data(cmptProduct(E-S));

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
    {
        D(i,j,k) = data[c++];
    }
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

    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_);
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_);

    List<floatScalar> data(cmptProduct(E-S)*n);

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

    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_);
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_);

    List<floatScalar> data(cmptProduct(E-S)*m*n);

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