#ifndef pMPI_H
#define pMPI_H

#include "label.H"
#include "scalar.H"
#include "VectorSpace.H"
#include "CellSpace.H"

#include "PstreamGlobals.H"

// Class that provides an MPI data type. When included, MPI must be included and
// linked.

namespace Foam
{

// Base class for templated types

template<class Type, template<class> class TType>
class pMPI<TType<Type>>
{
public:

    MPI_Datatype type;

    inline pMPI()
    {
        MPI_Type_contiguous
        (
            TType<Type>::nComponents*pTraits<Type>::nComponents,
            MPI_DOUBLE,
            &type
        );
        MPI_Type_commit(&type);
    }
};

// Non-templated primitive types

template<class Type>
class pMPI
{
public:

    MPI_Datatype type;

    inline pMPI()
    {
        MPI_Type_contiguous(pTraits<Type>::nComponents, MPI_DOUBLE, &type);
        MPI_Type_commit(&type);
    }
};

// Label specialization

template<>
class pMPI<label>
{
public:

    MPI_Datatype type;

    inline pMPI()
    {
        MPI_Type_contiguous(1, MPI_INT, &type);
        MPI_Type_commit(&type);
    }
};

}

#endif
