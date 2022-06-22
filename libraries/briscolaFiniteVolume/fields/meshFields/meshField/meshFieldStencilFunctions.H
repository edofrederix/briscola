#ifndef meshFieldStencilFunctions_H
#define meshFieldStencilFunctions_H

// Stencil products are defined by the '%' operator. They behave like matrix
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
    meshField<typename stencilProduct<Type1,Type2>::type,MeshType>&,
    const meshField<Type1,MeshType>&,
    const meshField<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshField<Type1,MeshType>&,
    const meshField<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshField<Type1,MeshType>>&,
    const meshField<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const meshField<Type1,MeshType>&,
    const tmp<meshField<Type2,MeshType>>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
Amul
(
    const tmp<meshField<Type1,MeshType>>&,
    const tmp<meshField<Type2,MeshType>>&
);

template<class Type, class Form, class MeshType>
void Amul
(
    meshField<typename stencilProduct<Type,Form>::type,MeshType>&,
    const meshField<Type,MeshType>&,
    const Form&
);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const meshField<Type,MeshType>&, const Form&);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
Amul(const tmp<meshField<Type,MeshType>>&, const Form&);

template<class Type, class Form, class MeshType>
void Amul
(
    meshField<typename stencilProduct<Form,Type>::type,MeshType>&,
    const Form&,
    const meshField<Type,MeshType>&
);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form&, const meshField<Type,MeshType>&);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
Amul(const Form&, const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
void rowSum
(
    meshField<scalar,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
rowSum(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
rowSum(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
scalarList matrixSum(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
scalarList matrixSum(const tmp<meshField<Type,MeshType>>&);

template<class Type, class MeshType>
void neighborSum
(
    meshField<scalar,MeshType>&,
    const meshField<Type,MeshType>&
);

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
neighborSum(const meshField<Type,MeshType>&);

template<class Type, class MeshType>
tmp<meshField<scalar,MeshType>>
neighborSum(const tmp<meshField<Type,MeshType>>&);

}

}

}

#endif
