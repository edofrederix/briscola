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
    List<autoPtr<std::ofstream>>& filePtrs,
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
        forAll(filePtrs, proc)
        {
            std::ofstream& file = filePtrs[proc]();

            const label size = (partitioned_ ? sizes[proc] : sum(sizes));

            if (header)
            {
                file<< name << " " << nComponents << " "
                    << size/nComponents << " float" << std::endl;
            }

            List<floatScalar> buffer;

            forAll(sizes, p)
            if (!partitioned_ || p == proc)
            {
                if (p == Pstream::masterNo())
                {
                    buffer.clear();
                    buffer = data;
                }
                else
                {
                    buffer.clear();
                    buffer.setSize(sizes[p]);

                    UIPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        p,
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
    List<autoPtr<std::ofstream>>& filePtrs,
    const meshDirection<Type,MeshType>& D
) const
{
    label c = 0;

    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_);
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_);
    const labelVector N = E - S;

    List<floatScalar> data(nStructured(N));

    for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (structured(i-S.x(),j-S.y(),k-S.z(),N))
                    data[c++] = D(i,j,k);

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::writeList
    (
        filePtrs,
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
    List<autoPtr<std::ofstream>>& filePtrs,
    const meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nComponents);

    label c = 0;

    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_);
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_);
    const labelVector N = E - S;

    List<floatScalar> data(nStructured(N)*n);

    for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (structured(i-S.x(),j-S.y(),k-S.z(),N))
                    for (label ii = 0; ii < n; ii++)
                        data[c++] = D(i,j,k)[ii];

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::writeList
    (
        filePtrs,
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
    List<autoPtr<std::ofstream>>& filePtrs,
    const meshDirection<Type,MeshType>& D
) const
{
    const label n(Type::nComponents);
    const label m(pTraits<typename Type::cmpt>::nComponents);

    label c = 0;

    const labelVector S = D.I().lower()-unitXYZ*label(ghosts_);
    const labelVector E = D.I().upper()+unitXYZ*label(ghosts_);
    const labelVector N = E - S;

    List<floatScalar> data(nStructured(N)*m*n);

    for (int i = S.x(); i < E.x(); i++)
        for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (structured(i-S.x(),j-S.y(),k-S.z(),N))
                    for (label ii = 0; ii < n; ii++)
                        for (label jj = 0; jj < m; jj++)
                            data[c++] = D(i,j,k)[ii][jj];

    const label tag =
        D.levelNum()*MeshType::numberOfDirections + D.directionNum();

    IO::writeList
    (
        filePtrs,
        data,
        D.mshLevel().mshField().name(),
        m*n,
        fvMsh_.time().writeFormat() == IOstream::ASCII,
        tag
    );
}

// Instantiate

#define WRITETYPEFIELD(FUNC,TYPE,MESHTYPE)                                     \
                                                                               \
template<>                                                                     \
void IO::writeField                                                            \
(                                                                              \
    List<autoPtr<std::ofstream>>& filePtrs,                                    \
    const meshDirection<TYPE,MESHTYPE>& D                                      \
) const                                                                        \
{                                                                              \
    FUNC(filePtrs,D);                                                          \
}                                                                              \
                                                                               \
template void IO::FUNC                                                         \
(                                                                              \
    List<autoPtr<std::ofstream>>&,                                             \
    const meshDirection<TYPE,MESHTYPE>&                                        \
) const;

WRITETYPEFIELD(writeScalarField,scalar,colocated)
WRITETYPEFIELD(writeScalarField,label,colocated)
WRITETYPEFIELD(writeArrayField,vector,colocated)
WRITETYPEFIELD(writeArrayField,tensor,colocated)
WRITETYPEFIELD(writeArrayField,diagTensor,colocated)
WRITETYPEFIELD(writeArrayField,symmTensor,colocated)
WRITETYPEFIELD(writeArrayField,sphericalTensor,colocated)
WRITETYPEFIELD(writeArrayField,faceScalar,colocated)
WRITETYPEFIELD(writeArrayField,edgeScalar,colocated)
WRITETYPEFIELD(writeArrayField,vertexScalar,colocated)
WRITETYPEFIELD(writeArrayArrayField,faceVector,colocated)
WRITETYPEFIELD(writeArrayArrayField,edgeVector,colocated)
WRITETYPEFIELD(writeArrayArrayField,vertexVector,colocated)

WRITETYPEFIELD(writeScalarField,scalar,staggered)
WRITETYPEFIELD(writeScalarField,label,staggered)
WRITETYPEFIELD(writeArrayField,vector,staggered)
WRITETYPEFIELD(writeArrayField,tensor,staggered)
WRITETYPEFIELD(writeArrayField,diagTensor,staggered)
WRITETYPEFIELD(writeArrayField,symmTensor,staggered)
WRITETYPEFIELD(writeArrayField,sphericalTensor,staggered)
WRITETYPEFIELD(writeArrayField,faceScalar,staggered)
WRITETYPEFIELD(writeArrayField,edgeScalar,staggered)
WRITETYPEFIELD(writeArrayField,vertexScalar,staggered)
WRITETYPEFIELD(writeArrayArrayField,faceVector,staggered)
WRITETYPEFIELD(writeArrayArrayField,edgeVector,staggered)
WRITETYPEFIELD(writeArrayArrayField,vertexVector,staggered)

}

}

}