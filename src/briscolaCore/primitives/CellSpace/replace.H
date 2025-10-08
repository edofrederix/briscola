#ifndef replace_H
#define replace_H

namespace Foam
{

// A special component replacement function is needed for CellSpace types. The
// default implementation below works for standard OpenFOAM primitives. It needs
// to be specialized for CellSpace types.

template<class Type, class CmptType>
inline void replace
(
    Type& p,
    const direction dir,
    const CmptType& v
)
{
    setComponent(p,dir) = v;
}

}

#endif
