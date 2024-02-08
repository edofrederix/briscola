#ifndef meshFieldStencilFunctions_H
#define meshFieldStencilFunctions_H

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type1, class Type2, class MeshType>
void rowProduct
(
    meshField<typename stencilProduct<Type1,Type2>::type,MeshType>&,
    const meshField<Type1,MeshType>&,
    const meshField<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const meshField<Type1,MeshType>&,
    const meshField<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshField<Type1,MeshType>>&,
    const meshField<Type2,MeshType>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const meshField<Type1,MeshType>&,
    const tmp<meshField<Type2,MeshType>>&
);

template<class Type1, class Type2, class MeshType>
tmp<meshField<typename stencilProduct<Type1,Type2>::type,MeshType>>
rowProduct
(
    const tmp<meshField<Type1,MeshType>>&,
    const tmp<meshField<Type2,MeshType>>&
);

template<class Type, class Form, class MeshType>
void rowProduct
(
    meshField<typename stencilProduct<Type,Form>::type,MeshType>&,
    const meshField<Type,MeshType>&,
    const Form&
);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const meshField<Type,MeshType>&, const Form&);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Type,Form>::type,MeshType>>
rowProduct(const tmp<meshField<Type,MeshType>>&, const Form&);

template<class Type, class Form, class MeshType>
void rowProduct
(
    meshField<typename stencilProduct<Form,Type>::type,MeshType>&,
    const Form&,
    const meshField<Type,MeshType>&
);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form&, const meshField<Type,MeshType>&);

template<class Type, class Form, class MeshType>
tmp<meshField<typename stencilProduct<Form,Type>::type,MeshType>>
rowProduct(const Form&, const tmp<meshField<Type,MeshType>>&);

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

}

}

}

#endif
