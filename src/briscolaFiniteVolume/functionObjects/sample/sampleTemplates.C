#include "sample.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

template<class Type>
void sample::appendScalarData
(
    const meshField<Type,colocated>& field,
    PtrList<scalarList>& data,
    wordList& headers
)
{
    pointInterpolator<colocated>& interp = interpPtr_();

    headers.append(field.name());
    data.append(new scalarList(move(interp(field))));
}

template<class Type>
void sample::appendArrayData
(
    const meshField<Type,colocated>& field,
    PtrList<scalarList>& data,
    wordList& headers
)
{
    pointInterpolator<colocated>& interp = interpPtr_();

    List<Type> interpData(move(interp(field)));

    const label n(Type::nComponents);

    for (int i = 0; i < n; i++)
    {
        data.append(new scalarList(interpData.size()));

        scalarList& datai = data[data.size()-1];

        forAll(datai, j)
        {
            datai[j] = interpData[j][i];
        }

        headers.append(field.name()+"."+pTraits<Type>::componentNames[i]);
    }
}

template<class Type>
void sample::appendArrayArrayData
(
    const meshField<Type,colocated>& field,
    PtrList<scalarList>& data,
    wordList& headers
)
{
    pointInterpolator<colocated>& interp = interpPtr_();

    List<Type> interpData(move(interp(field)));

    const label n(Type::nComponents);
    const label m(pTraits<typename Type::cmpt>::nComponents);

    for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
    {
        data.append(new scalarList(interpData.size()));

        scalarList& datai = data[data.size()-1];

        forAll(datai, k)
        {
            datai[k] = interpData[k][i][j];
        }

        headers.append
        (
            field.name()
          + "."
          + pTraits<Type>::componentNames[i]
          + "."
          + pTraits<typename Type::cmpt>::componentNames[j]
        );
    }
}

#define APPENDDATA(FUNC,TYPE)                                               \
                                                                            \
template<>                                                                  \
void sample::appendData                                                     \
(                                                                           \
    const meshField<TYPE,colocated>& field,                                 \
    PtrList<scalarList>& data,                                              \
    wordList& headers                                                       \
)                                                                           \
{                                                                           \
    this->FUNC(field,data,headers);                                         \
}                                                                           \
                                                                            \
template void sample::FUNC                                                  \
(                                                                           \
    const meshField<TYPE,colocated>& field,                                 \
    PtrList<scalarList>& data,                                              \
    wordList& headers                                                       \
);

APPENDDATA(appendScalarData,scalar)
APPENDDATA(appendArrayData,vector)
APPENDDATA(appendArrayData,tensor)
APPENDDATA(appendArrayData,sphericalTensor)
APPENDDATA(appendArrayData,symmTensor)
APPENDDATA(appendArrayData,diagTensor)
APPENDDATA(appendArrayData,faceScalar)
APPENDDATA(appendArrayArrayData,faceVector)

}

}

}

}
