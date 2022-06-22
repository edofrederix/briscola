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

// Default function for VectorSpace primitives

template<class Type, class MeshType>
void IO::readField
(
    vtkStructuredGrid * grid,
    meshField<Type,MeshType>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<Type,MeshType>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        field->GetTuple(c++, D(i,j,k).v_);
    }
}

// Instantiate

template void IO::readField
(
    vtkStructuredGrid*,
    colocatedVectorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    colocatedTensorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    colocatedSymmTensorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    colocatedSphericalTensorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    colocatedDiagTensorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    staggeredVectorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    staggeredTensorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    staggeredSymmTensorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    staggeredSphericalTensorField&,
    const label,
    const label
) const;

template void IO::readField
(
    vtkStructuredGrid*,
    staggeredDiagTensorField&,
    const label,
    const label
) const;

// Specialization for scalar type

template<>
void IO::readField
(
    vtkStructuredGrid * grid,
    meshField<scalar,colocated>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<scalar,colocated>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        D(i,j,k) = field->GetValue(c++);
    }
}

template<>
void IO::readField
(
    vtkStructuredGrid * grid,
    meshField<scalar,staggered>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<scalar,staggered>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        D(i,j,k) = field->GetValue(c++);
    }
}

// Specialization for label type

template<>
void IO::readField
(
    vtkStructuredGrid * grid,
    meshField<label,colocated>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<label,colocated>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        D(i,j,k) = field->GetValue(c++);
    }
}

template<>
void IO::readField
(
    vtkStructuredGrid * grid,
    meshField<label,staggered>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<label,staggered>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        D(i,j,k) = field->GetValue(c++);
    }
}

// Specialization for hex types

template<>
void IO::readField<hexScalar,colocated>
(
    vtkStructuredGrid * grid,
    meshField<hexScalar,colocated>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<hexScalar,colocated>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[6];

        field->GetTuple(c++, v);

        D(i,j,k).left() = v[0];
        D(i,j,k).right() = v[1];
        D(i,j,k).bottom() = v[2];
        D(i,j,k).top() = v[3];
        D(i,j,k).aft() = v[4];
        D(i,j,k).fore() = v[5];
    }
}

template<>
void IO::readField<hexScalar,staggered>
(
    vtkStructuredGrid * grid,
    meshField<hexScalar,staggered>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<hexScalar,staggered>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[6];

        field->GetTuple(c++, v);

        D(i,j,k).left() = v[0];
        D(i,j,k).right() = v[1];
        D(i,j,k).bottom() = v[2];
        D(i,j,k).top() = v[3];
        D(i,j,k).aft() = v[4];
        D(i,j,k).fore() = v[5];
    }
}

template<>
void IO::readField<hexVector,colocated>
(
    vtkStructuredGrid * grid,
    meshField<hexVector,colocated>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<hexVector,colocated>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[18];

        field->GetTuple(c++, v);

        D(i,j,k).left().x() = v[0];
        D(i,j,k).left().y() = v[1];
        D(i,j,k).left().z() = v[2];

        D(i,j,k).right().x() = v[3];
        D(i,j,k).right().y() = v[4];
        D(i,j,k).right().z() = v[5];

        D(i,j,k).bottom().x() = v[6];
        D(i,j,k).bottom().y() = v[7];
        D(i,j,k).bottom().z() = v[8];

        D(i,j,k).top().x() = v[9];
        D(i,j,k).top().y() = v[10];
        D(i,j,k).top().z() = v[11];

        D(i,j,k).aft().x() = v[12];
        D(i,j,k).aft().y() = v[13];
        D(i,j,k).aft().z() = v[14];

        D(i,j,k).fore().x() = v[15];
        D(i,j,k).fore().y() = v[16];
        D(i,j,k).fore().z() = v[17];
    }
}

template<>
void IO::readField<hexVector,staggered>
(
    vtkStructuredGrid * grid,
    meshField<hexVector,staggered>& f,
    const label l,
    const label d
) const
{
    const word name(f.name());

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

    meshDirection<hexVector,staggered>& D = f[l][d];

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[18];

        field->GetTuple(c++, v);

        D(i,j,k).left().x() = v[0];
        D(i,j,k).left().y() = v[1];
        D(i,j,k).left().z() = v[2];

        D(i,j,k).right().x() = v[3];
        D(i,j,k).right().y() = v[4];
        D(i,j,k).right().z() = v[5];

        D(i,j,k).bottom().x() = v[6];
        D(i,j,k).bottom().y() = v[7];
        D(i,j,k).bottom().z() = v[8];

        D(i,j,k).top().x() = v[9];
        D(i,j,k).top().y() = v[10];
        D(i,j,k).top().z() = v[11];

        D(i,j,k).aft().x() = v[12];
        D(i,j,k).aft().y() = v[13];
        D(i,j,k).aft().z() = v[14];

        D(i,j,k).fore().x() = v[15];
        D(i,j,k).fore().y() = v[16];
        D(i,j,k).fore().z() = v[17];
    }
}

}

}

}