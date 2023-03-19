#ifndef faceLabel_H
#define faceLabel_H

#include "label.H"
#include "FaceSpace.H"
#include "contiguous.H"

namespace Foam
{

typedef FaceSpace<label> faceLabel;

template<>
inline bool contiguous<faceLabel>() { return true; }

// Specialization of CellSpaceCast for stencils

template<>
class CellSpaceCast<faceLabel, stencil, label>
{
public:
    typedef stencil type;
};

template<>
class CellSpaceCast<stencil, faceLabel, label>
{
public:
    typedef stencil type;
};

}

#endif