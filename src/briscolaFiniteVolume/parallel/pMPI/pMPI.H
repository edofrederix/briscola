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

    ~pMPI()
    {
        MPI_Type_free(&type);
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

    ~pMPI()
    {
        MPI_Type_free(&type);
    }
};

// Cell space specializations

template<class Type>
class csMPI
{
public:

    MPI_Datatype type;

    inline csMPI()
    {
        MPI_Type_contiguous
        (
            Type::nCsComponents*Type::nComponents,
            MPI_DOUBLE,
            &type
        );
        MPI_Type_commit(&type);
    }

    ~csMPI()
    {
        MPI_Type_free(&type);
    }
};

// Face space

template<> class pMPI<faceScalar> : public csMPI<faceScalar::csType> {};
template<> class pMPI<faceVector> : public csMPI<faceVector::csType> {};

// Edge space

template<> class pMPI<edgeScalar> : public csMPI<edgeScalar::csType> {};
template<> class pMPI<edgeVector> : public csMPI<edgeVector::csType> {};

// Vertex space

template<> class pMPI<vertexScalar> : public csMPI<vertexScalar::csType> {};
template<> class pMPI<vertexVector> : public csMPI<vertexVector::csType> {};

// Stencil space

template<> class pMPI<stencil> : public csMPI<stencil::csType> {};
template<> class pMPI<diagStencil> : public csMPI<diagStencil::csType> {};

}

#endif
