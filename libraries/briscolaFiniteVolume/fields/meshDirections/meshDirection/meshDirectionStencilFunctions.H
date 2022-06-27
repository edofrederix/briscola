#ifndef meshDirectionStencilFunctions_H
#define meshDirectionStencilFunctions_H

#include "meshDirectionStencilProduct.H"
#include "meshDirectionDiagStencilProduct.H"

// Stencil products are defined by the Amul function. They behave like matrix
// multiplications, but require non-local operations. Those operations are
// defined in the Amul functions operating on a mesh level

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void Amul
(
    meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>&,
    const meshDirection<Type1,MeshType>&,
    const meshDirection<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshDirection<Type1,MeshType>&,
    const meshDirection<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshDirection<Type1,MeshType>>&,
    const meshDirection<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshDirection<Type1,MeshType>&,
    const tmp<meshDirection<Type2,MeshType>>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshDirection<Type1,MeshType>>&,
    const tmp<meshDirection<Type2,MeshType>>&
);

template<class Type, class Form, class MeshType>
void Amul
(
    meshDirection<typename stencilProduct<Type,Form>::type,MeshType>&,
    const meshDirection<Type,MeshType>&,
    const Form&
);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const meshDirection<Type,MeshType>&, const Form&);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const tmp<meshDirection<Type,MeshType>>&, const Form&);

template<class Type, class Form, class MeshType>
void Amul
(
    meshDirection<typename stencilProduct<Form,Type>::type,MeshType>&,
    const Form&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form&, const meshDirection<Type,MeshType>&);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form&, const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
void rowSum
(
    meshDirection<scalar,MeshType>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
rowSum(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
rowSum(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
scalar matrixSum(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
scalar matrixSum(const tmp<meshDirection<Type,MeshType>>&);

template<class Type, class MeshType>
void neighborSum
(
    meshDirection<scalar,MeshType>&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
neighborSum(const meshDirection<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshDirection<scalar,MeshType>>
neighborSum(const tmp<meshDirection<Type,MeshType>>&);

}

}

}

#endif
