#ifndef faceScalar_H
#define faceScalar_H

#include "scalar.H"
#include "FaceSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef FaceSpace<scalar> faceScalar;

template<>
inline bool contiguous<faceScalar>() { return true; }

// Specialization of CellSpaceCast for stencils

template<>
class CellSpaceCast<faceScalar, stencil, scalar>
{
public:
    typedef stencil type;
};

template<>
class CellSpaceCast<stencil, faceScalar, scalar>
{
public:
    typedef stencil type;
};

}

#include "faceScalarI.H"

#endif