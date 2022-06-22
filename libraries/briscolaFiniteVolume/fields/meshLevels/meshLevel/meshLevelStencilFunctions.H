#ifndef meshLevelStencilFunctions_H
#define meshLevelStencilFunctions_H

#include "meshDirectionStencilFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void Amul
(
    meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>&,
    const meshLevel<Type1,MeshType>&,
    const meshLevel<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshLevel<Type1,MeshType>&,
    const meshLevel<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshLevel<Type1,MeshType>>&,
    const meshLevel<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshLevel<Type1,MeshType>&,
    const tmp<meshLevel<Type2,MeshType>>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshLevel<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshLevel<Type1,MeshType>>&,
    const tmp<meshLevel<Type2,MeshType>>&
);

template<class Type, class Form, class MeshType>
void Amul
(
    meshLevel<typename stencilProduct<Type,Form>::type,MeshType>&,
    const meshLevel<Type,MeshType>&,
    const Form&
);

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const meshLevel<Type,MeshType>&, const Form&);

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const tmp<meshLevel<Type,MeshType>>&, const Form&);

template<class Type, class Form, class MeshType>
void Amul
(
    meshLevel<typename stencilProduct<Form,Type>::type,MeshType>&,
    const Form&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form&, const meshLevel<Type,MeshType>&);

template<class Type, class Form, class MeshType>
tmp<meshLevel<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form&, const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void rowSum
(
    meshLevel<scalar,MeshType>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
rowSum(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
rowSum(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
scalarList matrixSum(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
scalarList matrixSum(const tmp<meshLevel<Type,MeshType>>&);

template<class Type, class MeshType>
void neighborSum
(
    meshLevel<scalar,MeshType>&,
    const meshLevel<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
neighborSum(const meshLevel<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshLevel<scalar,MeshType>>
neighborSum(const tmp<meshLevel<Type,MeshType>>&);

}

}

}

#endif
