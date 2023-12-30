#ifndef pMPI_H
#define pMPI_H

#include "label.H"
#include "scalar.H"
#include "VectorSpace.H"
#include "CellSpace.H"

#include <mpi.h>

// Class that provides an MPI data type. When included, MPI must be included and
// linked. By default, the MPI data type is generated from nComponents and
// MPI_DOUBLE.

namespace Foam
{

// Base class

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

// Scalar specialization

template<>
class pMPI<scalar>
{
public:

    MPI_Datatype type;

    inline pMPI()
    {
        MPI_Type_contiguous(1, MPI_DOUBLE, &type);
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
