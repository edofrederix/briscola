#include "IO.H"
#include "colocatedFields.H"
#include "staggeredFields.H"

#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "vtkCellData.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void IO::readScalarField
(
    vtkStructuredGrid * grid,
    meshDirection<Type,MeshType>& D
) const
{
    const word name(D.mshLevel().mshField().name());

    if (!grid->GetCellData()->HasArray(name.c_str()))
    {
        FatalErrorInFunction
            << "Could not find VTK data for field " << name << endl
            << exit(FatalError);
    }

    vtkSmartPointer<vtkIntArray> field =
        vtkIntArray::FastDownCast
        (
            grid->GetCellData()->GetAbstractArray(name.c_str())
        );

    label c = 0;

    for (int k = 0; k < D.B().n()-2; k++)
    for (int j = 0; j < D.B().m()-2; j++)
    for (int i = 0; i < D.B().l()-2; i++)
    {
        D(i,j,k) = field->GetValue(c++);
    }
}

template<class Type, class MeshType>
void IO::readArrayField
(
    vtkStructuredGrid * grid,
    meshDirection<Type,MeshType>& D
) const
{
    const word name(D.mshLevel().mshField().name());

    if (!grid->GetCellData()->HasArray(name.c_str()))
    {
        FatalErrorInFunction
            << "Could not find VTK data for field " << name << endl
            << exit(FatalError);
    }

    vtkSmartPointer<vtkDoubleArray> field =
        vtkDoubleArray::FastDownCast
        (
            grid->GetCellData()->GetAbstractArray(name.c_str())
        );

    label c = 0;

    for (int k = 0; k < D.B().n()-2; k++)
    for (int j = 0; j < D.B().m()-2; j++)
    for (int i = 0; i < D.B().l()-2; i++)
    {
        field->GetTuple(c++, D(i,j,k).v_);
    }
}

template<class Type, class MeshType>
void IO::readArrayArrayField
(
    vtkStructuredGrid * grid,
    meshDirection<Type,MeshType>& D
) const
{
    const word name(D.mshLevel().mshField().name());

    const label m(Type::nComponents);
    const label n(pTraits<typename Type::cmpt>::nComponents);

    if (!grid->GetCellData()->HasArray(name.c_str()))
    {
        FatalErrorInFunction
            << "Could not find VTK data for field " << name << endl
            << exit(FatalError);
    }

    vtkSmartPointer<vtkDoubleArray> field =
        vtkDoubleArray::FastDownCast
        (
            grid->GetCellData()->GetAbstractArray(name.c_str())
        );

    label c = 0;

    for (int k = 0; k < D.B().n()-2; k++)
    for (int j = 0; j < D.B().m()-2; j++)
    for (int i = 0; i < D.B().l()-2; i++)
    {
        typename Type::cmptCmpt ar[m*n];
        field->GetTuple(c++, ar);

        for (label ii = 0; ii < m; ii++)
        {
            for (label jj = 0; jj < n; jj++)
            {
                D(i,j,k)[ii][jj] = ar[ii*n+jj];
            }
        }
    }
}

// Instantiate

#define READTYPEFIELD(FUNC,TYPE,MESHTYPE)                                       \
                                                                                \
template<>                                                                      \
void IO::readField                                                              \
(                                                                               \
    vtkStructuredGrid * grid,                                                   \
    meshDirection<TYPE,MESHTYPE>& D                                             \
) const                                                                         \
{                                                                               \
    FUNC(grid,D);                                                               \
}                                                                               \
                                                                                \
template void IO::FUNC                                                          \
(                                                                               \
    vtkStructuredGrid*,                                                         \
    meshDirection<TYPE,MESHTYPE>&                                               \
) const;

READTYPEFIELD(readScalarField,scalar,colocated)
READTYPEFIELD(readScalarField,label,colocated)
READTYPEFIELD(readArrayField,vector,colocated)
READTYPEFIELD(readArrayField,tensor,colocated)
READTYPEFIELD(readArrayField,diagTensor,colocated)
READTYPEFIELD(readArrayField,symmTensor,colocated)
READTYPEFIELD(readArrayField,sphericalTensor,colocated)
READTYPEFIELD(readArrayField,faceScalar,colocated)
READTYPEFIELD(readArrayArrayField,faceVector,colocated)

READTYPEFIELD(readScalarField,scalar,staggered)
READTYPEFIELD(readScalarField,label,staggered)
READTYPEFIELD(readArrayField,vector,staggered)
READTYPEFIELD(readArrayField,tensor,staggered)
READTYPEFIELD(readArrayField,diagTensor,staggered)
READTYPEFIELD(readArrayField,symmTensor,staggered)
READTYPEFIELD(readArrayField,sphericalTensor,staggered)
READTYPEFIELD(readArrayField,faceScalar,staggered)
READTYPEFIELD(readArrayArrayField,faceVector,staggered)

#undef READTYPEFIELD

}

}

}