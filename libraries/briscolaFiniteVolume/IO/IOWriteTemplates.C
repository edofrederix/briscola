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
void IO::writeField
(
    vtkStructuredGrid * grid,
    const meshField<Type,MeshType>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const meshDirection<Type,MeshType>& D = f[l][d];

    field->SetName(f.name().c_str());

    const int n = pTraits<Type>::nComponents;

    field->SetNumberOfComponents(n);
    field->SetNumberOfTuples(D.size());

    for (label c = 0; c < n; c++)
    {
        field->SetComponentName
        (
            c,
            pTraits<Type>::componentNames[c]
        );
    }

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        field->SetTuple(c++, D(i,j,k).v_);
    }

    grid->GetCellData()->AddArray(field);
}

// Instantiate

template void IO::writeField
(
    vtkStructuredGrid*,
    const colocatedVectorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const colocatedTensorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const colocatedSymmTensorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const colocatedSphericalTensorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const colocatedDiagTensorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const staggeredVectorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const staggeredTensorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const staggeredSymmTensorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const staggeredSphericalTensorField&,
    const label,
    const label
) const;

template void IO::writeField
(
    vtkStructuredGrid*,
    const staggeredDiagTensorField&,
    const label,
    const label
) const;

// Specialization for scalar type

template<>
void IO::writeField
(
    vtkStructuredGrid * grid,
    const meshField<scalar,colocated>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const meshDirection<scalar,colocated>& D = f[l][d];

    field->SetName(f.name().c_str());
    field->SetNumberOfValues(D.size());

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        field->SetValue(c++, D(i,j,k));
    }

    grid->GetCellData()->AddArray(field);
}

template<>
void IO::writeField
(
    vtkStructuredGrid * grid,
    const meshField<scalar,staggered>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const meshDirection<scalar,staggered>& D = f[l][d];

    field->SetName(f.name().c_str());
    field->SetNumberOfValues(D.size());

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        field->SetValue(c++, D(i,j,k));
    }

    grid->GetCellData()->AddArray(field);
}

// Specialization for label type

template<>
void IO::writeField
(
    vtkStructuredGrid * grid,
    const meshField<label,colocated>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkIntArray> field =
        vtkSmartPointer<vtkIntArray>::New();

    const meshDirection<label,colocated>& D = f[l][d];

    field->SetName(f.name().c_str());
    field->SetNumberOfValues(D.size());

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        field->SetValue(c++, D(i,j,k));
    }

    grid->GetCellData()->AddArray(field);
}

template<>
void IO::writeField
(
    vtkStructuredGrid * grid,
    const meshField<label,staggered>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkIntArray> field =
        vtkSmartPointer<vtkIntArray>::New();

    const meshDirection<label,staggered>& D = f[l][d];

    field->SetName(f.name().c_str());
    field->SetNumberOfValues(D.size());

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        field->SetValue(c++, D(i,j,k));
    }

    grid->GetCellData()->AddArray(field);
}

// Specialization for hex types

template<>
void IO::writeField<hexScalar,colocated>
(
    vtkStructuredGrid * grid,
    const meshField<hexScalar,colocated>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const meshDirection<hexScalar,colocated>& D = f[l][d];

    field->SetName(f.name().c_str());

    field->SetNumberOfComponents(6);
    field->SetNumberOfTuples(D.size());

    field->SetComponentName(0, "l");
    field->SetComponentName(1, "r");
    field->SetComponentName(2, "b");
    field->SetComponentName(3, "t");
    field->SetComponentName(4, "a");
    field->SetComponentName(5, "f");

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[6];

        v[0] = D(i,j,k).left();
        v[1] = D(i,j,k).right();
        v[2] = D(i,j,k).bottom();
        v[3] = D(i,j,k).top();
        v[4] = D(i,j,k).aft();
        v[5] = D(i,j,k).fore();

        field->SetTuple(c++, v);
    }

    grid->GetCellData()->AddArray(field);
}

template<>
void IO::writeField<hexScalar,staggered>
(
    vtkStructuredGrid * grid,
    const meshField<hexScalar,staggered>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const meshDirection<hexScalar,staggered>& D = f[l][d];

    field->SetName(f.name().c_str());

    field->SetNumberOfComponents(6);
    field->SetNumberOfTuples(D.size());

    field->SetComponentName(0, "l");
    field->SetComponentName(1, "r");
    field->SetComponentName(2, "b");
    field->SetComponentName(3, "t");
    field->SetComponentName(4, "a");
    field->SetComponentName(5, "f");

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[6];

        v[0] = D(i,j,k).left();
        v[1] = D(i,j,k).right();
        v[2] = D(i,j,k).bottom();
        v[3] = D(i,j,k).top();
        v[4] = D(i,j,k).aft();
        v[5] = D(i,j,k).fore();

        field->SetTuple(c++, v);
    }

    grid->GetCellData()->AddArray(field);
}

template<>
void IO::writeField<hexVector,colocated>
(
    vtkStructuredGrid * grid,
    const meshField<hexVector,colocated>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const meshDirection<hexVector,colocated>& D = f[l][d];

    field->SetName(f.name().c_str());

    field->SetNumberOfComponents(18);
    field->SetNumberOfTuples(D.size());

    field->SetComponentName(0, "lx");
    field->SetComponentName(1, "ly");
    field->SetComponentName(2, "lz");

    field->SetComponentName(3, "rx");
    field->SetComponentName(4, "ry");
    field->SetComponentName(5, "rz");

    field->SetComponentName(6, "bx");
    field->SetComponentName(7, "by");
    field->SetComponentName(8, "bz");

    field->SetComponentName(9, "tx");
    field->SetComponentName(10, "ty");
    field->SetComponentName(11, "tz");

    field->SetComponentName(12, "ax");
    field->SetComponentName(13, "ay");
    field->SetComponentName(14, "az");

    field->SetComponentName(15, "fx");
    field->SetComponentName(16, "fy");
    field->SetComponentName(17, "fz");

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[18];

        v[0] = D(i,j,k).left().x();
        v[1] = D(i,j,k).left().y();
        v[2] = D(i,j,k).left().z();

        v[3] = D(i,j,k).right().x();
        v[4] = D(i,j,k).right().y();
        v[5] = D(i,j,k).right().z();

        v[6] = D(i,j,k).bottom().x();
        v[7] = D(i,j,k).bottom().y();
        v[8] = D(i,j,k).bottom().z();

        v[9] = D(i,j,k).top().x();
        v[10] = D(i,j,k).top().y();
        v[11] = D(i,j,k).top().z();

        v[12] = D(i,j,k).aft().x();
        v[13] = D(i,j,k).aft().y();
        v[14] = D(i,j,k).aft().z();

        v[15] = D(i,j,k).fore().x();
        v[16] = D(i,j,k).fore().y();
        v[17] = D(i,j,k).fore().z();

        field->SetTuple(c++, v);
    }

    grid->GetCellData()->AddArray(field);
}

template<>
void IO::writeField<hexVector,staggered>
(
    vtkStructuredGrid * grid,
    const meshField<hexVector,staggered>& f,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkDoubleArray> field =
        vtkSmartPointer<vtkDoubleArray>::New();

    const meshDirection<hexVector,staggered>& D = f[l][d];

    field->SetName(f.name().c_str());

    field->SetNumberOfComponents(18);
    field->SetNumberOfTuples(D.size());

    field->SetComponentName(0, "lx");
    field->SetComponentName(1, "ly");
    field->SetComponentName(2, "lz");

    field->SetComponentName(3, "rx");
    field->SetComponentName(4, "ry");
    field->SetComponentName(5, "rz");

    field->SetComponentName(6, "bx");
    field->SetComponentName(7, "by");
    field->SetComponentName(8, "bz");

    field->SetComponentName(9, "tx");
    field->SetComponentName(10, "ty");
    field->SetComponentName(11, "tz");

    field->SetComponentName(12, "ax");
    field->SetComponentName(13, "ay");
    field->SetComponentName(14, "az");

    field->SetComponentName(15, "fx");
    field->SetComponentName(16, "fy");
    field->SetComponentName(17, "fz");

    label c = 0;

    forAllCellsReversed(D, i, j, k)
    {
        double v[18];

        v[0] = D(i,j,k).left().x();
        v[1] = D(i,j,k).left().y();
        v[2] = D(i,j,k).left().z();

        v[3] = D(i,j,k).right().x();
        v[4] = D(i,j,k).right().y();
        v[5] = D(i,j,k).right().z();

        v[6] = D(i,j,k).bottom().x();
        v[7] = D(i,j,k).bottom().y();
        v[8] = D(i,j,k).bottom().z();

        v[9] = D(i,j,k).top().x();
        v[10] = D(i,j,k).top().y();
        v[11] = D(i,j,k).top().z();

        v[12] = D(i,j,k).aft().x();
        v[13] = D(i,j,k).aft().y();
        v[14] = D(i,j,k).aft().z();

        v[15] = D(i,j,k).fore().x();
        v[16] = D(i,j,k).fore().y();
        v[17] = D(i,j,k).fore().z();

        field->SetTuple(c++, v);
    }

    grid->GetCellData()->AddArray(field);
}

}

}

}

