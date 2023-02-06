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

void IO::writeList
(
    autoPtr<std::ofstream>& filePtr,
    List<floatScalar>& data,
    const word name,
    const label nComponents,
    const bool ascii,
    const label tag,
    const bool header
) const
{
    labelList sizes(Pstream::nProcs());
    sizes[Pstream::myProcNo()] = data.size();
    Pstream::gatherList(sizes);

    if (Pstream::master())
    {
        std::ofstream& file = filePtr();

        if (header)
        {
            file<< name << " " << nComponents << " "
                << sum(sizes)/nComponents << " float" << std::endl;
        }

        List<floatScalar> buffer;

        forAll(sizes, proc)
        {
            if (proc == Pstream::masterNo())
            {
                buffer.clear();
                buffer = data;
            }
            else
            {
                buffer.clear();
                buffer.setSize(sizes[proc]);

                UIPstream::read
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    reinterpret_cast<char*>(buffer.begin()),
                    buffer.byteSize(),
                    tag,
                    UPstream::worldComm
                );
            }

            if (ascii)
            {
                forAll(buffer, i)
                {
                    file<< buffer[i] << " ";
                }
            }
            else
            {
                #ifdef LITTLEENDIAN
                swapWords
                (
                    buffer.byteSize()/sizeof(label),
                    reinterpret_cast<label*>(buffer.begin())
                );
                #endif

                file.write
                (
                    reinterpret_cast<char*>(buffer.begin()),
                    buffer.byteSize()
                );
            }
        }

        file<< std::endl;
    }
    else
    {
        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo(),
            reinterpret_cast<char*>(data.begin()),
            data.byteSize(),
            tag,
            UPstream::worldComm
        );
    }

    Pstream::waitRequests();
}

template<class Type, class MeshType>
void IO::writeScalarField
(
    autoPtr<std::ofstream>& filePtr,
    const meshDirection<Type,MeshType>& D
) const
{
    List<floatScalar> data(D.size());

    label c = 0;

    forAllCells(D, i, j, k)
        data[c++] = D(i,j,k);

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::writeList
    (
        filePtr,
        data,
        D.mshLevel().mshField().name(),
        1,
        fvMsh_.time().writeFormat() == IOstream::ASCII,
        tag
    );
}

template<class Type, class MeshType>
void IO::writeArrayField
(
    autoPtr<std::ofstream>& filePtr,
    const meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nComponents);

    List<floatScalar> data(D.size()*n);

    label c = 0;

    forAllCells(D, i, j, k)
        for (label ii = 0; ii < n; ii++)
            data[c++] = D(i,j,k)[ii];

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::writeList
    (
        filePtr,
        data,
        D.mshLevel().mshField().name(),
        n,
        fvMsh_.time().writeFormat() == IOstream::ASCII,
        tag
    );
}

template<class Type, class MeshType>
void IO::writeArrayArrayField
(
    autoPtr<std::ofstream>& filePtr,
    const meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nComponents);
    const label m(pTraits<typename Type::cmpt>::nComponents);

    List<floatScalar> data(D.size()*m*n);

    label c = 0;
    forAllCells(D, i, j, k)
        for (label ii = 0; ii < n; ii++)
            for (label jj = 0; jj < m; jj++)
                data[c++] = D(i,j,k)[ii][jj];

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::writeList
    (
        filePtr,
        data,
        D.mshLevel().mshField().name(),
        m*n,
        fvMsh_.time().writeFormat() == IOstream::ASCII,
        tag
    );
}

// Instantiate

#define WRITETYPEFIELD(FUNC,TYPE,MESHTYPE)                                      \
                                                                                \
template<>                                                                      \
void IO::writeField                                                             \
(                                                                               \
    autoPtr<std::ofstream>& filePtr,                                            \
    const meshDirection<TYPE,MESHTYPE>& D                                       \
) const                                                                         \
{                                                                               \
    FUNC(filePtr,D);                                                            \
}                                                                               \
                                                                                \
template void IO::FUNC                                                          \
(                                                                               \
    autoPtr<std::ofstream>&,                                                    \
    const meshDirection<TYPE,MESHTYPE>&                                         \
) const;

WRITETYPEFIELD(writeScalarField,scalar,colocated)
WRITETYPEFIELD(writeScalarField,label,colocated)
WRITETYPEFIELD(writeArrayField,vector,colocated)
WRITETYPEFIELD(writeArrayField,tensor,colocated)
WRITETYPEFIELD(writeArrayField,diagTensor,colocated)
WRITETYPEFIELD(writeArrayField,symmTensor,colocated)
WRITETYPEFIELD(writeArrayField,sphericalTensor,colocated)
WRITETYPEFIELD(writeArrayField,faceScalar,colocated)
WRITETYPEFIELD(writeArrayArrayField,faceVector,colocated)

WRITETYPEFIELD(writeScalarField,scalar,staggered)
WRITETYPEFIELD(writeScalarField,label,staggered)
WRITETYPEFIELD(writeArrayField,vector,staggered)
WRITETYPEFIELD(writeArrayField,tensor,staggered)
WRITETYPEFIELD(writeArrayField,diagTensor,staggered)
WRITETYPEFIELD(writeArrayField,symmTensor,staggered)
WRITETYPEFIELD(writeArrayField,sphericalTensor,staggered)
WRITETYPEFIELD(writeArrayField,faceScalar,staggered)
WRITETYPEFIELD(writeArrayArrayField,faceVector,staggered)

}

}

}