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
void IO::writeScalarField
(
    vtkStructuredGrid * grid,
    const meshDirection<Type,MeshType>& D
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    field->SetName(D.mshLevel().mshField().name().c_str());
    field->SetNumberOfValues(cmptProduct(D.B().N()-2*unitXYZ));

    // Also write inactive internal cells. Use the block shape and subtract
    // ghost cells.

    label c = 0;

    for (int k = 0; k < D.B().n()-2; k++)
    for (int j = 0; j < D.B().m()-2; j++)
    for (int i = 0; i < D.B().l()-2; i++)
    {
        field->SetValue(c++, D(i,j,k));
    }

    grid->GetCellData()->AddArray(field);
}

template<class Type, class MeshType>
void IO::writeArrayField
(
    vtkStructuredGrid * grid,
    const meshDirection<Type,MeshType>& D
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const label n(Type::nComponents);

    field->SetName(D.mshLevel().mshField().name().c_str());
    field->SetNumberOfComponents(n);
    field->SetNumberOfTuples(cmptProduct(D.B().N()-2*unitXYZ));

    for (label i = 0; i < n; i++)
    {
        field->SetComponentName
        (
            i,
            Type::componentNames[i]
        );
    }

    // Also write inactive internal cells. Use the block shape and subtract
    // ghost cells.

    label c = 0;

    for (int k = 0; k < D.B().n()-2; k++)
    for (int j = 0; j < D.B().m()-2; j++)
    for (int i = 0; i < D.B().l()-2; i++)
    {
        field->SetTuple(c++, D(i,j,k).v_);
    }

    grid->GetCellData()->AddArray(field);
}

template<class Type, class MeshType>
void IO::writeArrayArrayField
(
    vtkStructuredGrid * grid,
    const meshDirection<Type,MeshType>& D
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const label m(Type::nComponents);
    const label n(pTraits<typename Type::cmpt>::nComponents);

    field->SetName(D.mshLevel().mshField().name().c_str());
    field->SetNumberOfComponents(m*n);
    field->SetNumberOfTuples(cmptProduct(D.B().N()-2*unitXYZ));

    for (label i = 0; i < m; i++)
    {
        for (label j = 0; j < n; j++)
        {
            field->SetComponentName
            (
                i*n+j,
                (
                    word(Type::componentNames[i]) + "." +
                    word(pTraits<typename Type::cmpt>::componentNames[j])
                ).c_str()
            );
        }
    }

    // Also write inactive internal cells. Use the block shape and subtract
    // ghost cells.

    label c = 0;

    for (int k = 0; k < D.B().n()-2; k++)
    for (int j = 0; j < D.B().m()-2; j++)
    for (int i = 0; i < D.B().l()-2; i++)
    {
        typename Type::cmptCmpt ar[m*n];

        for (label ii = 0; ii < m; ii++)
        {
            for (label jj = 0; jj < n; jj++)
            {
                ar[ii*n+jj] = D(i,j,k)[ii][jj];
            }
        }

        field->SetTuple(c++, ar);
    }

    grid->GetCellData()->AddArray(field);
}

// Instantiate

#define WRITETYPEFIELD(FUNC,TYPE,MESHTYPE)                                      \
                                                                                \
template<>                                                                      \
void IO::writeField                                                             \
(                                                                               \
    vtkStructuredGrid * grid,                                                   \
    const meshDirection<TYPE,MESHTYPE>& D                                       \
) const                                                                         \
{                                                                               \
    FUNC(grid,D);                                                               \
}                                                                               \
                                                                                \
template void IO::FUNC                                                          \
(                                                                               \
    vtkStructuredGrid*,                                                         \
    const meshDirection<TYPE,MESHTYPE>&                                         \
) const;

WRITETYPEFIELD(writeScalarField,scalar,colocated)
WRITETYPEFIELD(writeScalarField,label,colocated)
WRITETYPEFIELD(writeArrayField,vector,colocated)
WRITETYPEFIELD(writeArrayField,tensor,colocated)
WRITETYPEFIELD(writeArrayField,diagTensor,colocated)
WRITETYPEFIELD(writeArrayField,symmTensor,colocated)
WRITETYPEFIELD(writeArrayField,sphericalTensor,colocated)
WRITETYPEFIELD(writeArrayField,faceScalar,colocated)
WRITETYPEFIELD(writeArrayArrayField,faceVector,colocated)

WRITETYPEFIELD(writeScalarField,scalar,staggered)
WRITETYPEFIELD(writeScalarField,label,staggered)
WRITETYPEFIELD(writeArrayField,vector,staggered)
WRITETYPEFIELD(writeArrayField,tensor,staggered)
WRITETYPEFIELD(writeArrayField,diagTensor,staggered)
WRITETYPEFIELD(writeArrayField,symmTensor,staggered)
WRITETYPEFIELD(writeArrayField,sphericalTensor,staggered)
WRITETYPEFIELD(writeArrayField,faceScalar,staggered)
WRITETYPEFIELD(writeArrayArrayField,faceVector,staggered)

}

}

}