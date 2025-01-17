#ifndef meshDirectionStencilFunctions_H
#define meshDirectionStencilFunctions_H

#include "meshDirectionStencilOperations.H"
#include "meshDirectionDiagStencilOperations.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void rowProduct
(
    meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>&,
    const meshDirection<Type1,MeshType>&,
    const meshDirection<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const meshDirection<Type1,MeshType>&,
    const meshDirection<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshDirection<Type1,MeshType>>&,
    const meshDirection<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const meshDirection<Type1,MeshType>&,
    const tmp<meshDirection<Type2,MeshType>>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshDirection<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshDirection<Type1,MeshType>>&,
    const tmp<meshDirection<Type2,MeshType>>&
);

template<class Type, class Form, class MeshType>
void rowProduct
(
    meshDirection<typename stencilProduct<Type,Form>::type,MeshType>&,
    const meshDirection<Type,MeshType>&,
    const Form&
);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const meshDirection<Type,MeshType>&, const Form&);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const tmp<meshDirection<Type,MeshType>>&, const Form&);

template<class Type, class Form, class MeshType>
void rowProduct
(
    meshDirection<typename stencilProduct<Form,Type>::type,MeshType>&,
    const Form&,
    const meshDirection<Type,MeshType>&
);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form&, const meshDirection<Type,MeshType>&);

template<class Type, class Form, class MeshType>
tmp<meshDirection<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form&, const tmp<meshDirection<Type,MeshType>>&);

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

}

}

}

#endif
