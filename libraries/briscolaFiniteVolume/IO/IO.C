#include "IO.H"
#include "OSspecific.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"
#include "SortableList.H"
#include "colocatedFields.H"
#include "staggeredFields.H"

#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkHexahedron.h"
#include "vtkXMLStructuredGridWriter.h"
#include "vtkXMLStructuredGridReader.h"
#include "vtkCellData.h"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(IO, 0);

word IO::dataFileNamePrefix("briscola");
word IO::timeDataFileName("timeData");
word IO::XMLStructuredGridExtension("vts");

template<class MeshType>
void IO::writeData(const word timeName, const label l)
{
    // Set points

    const partLevelPoints& p = fvMsh_[l].points();

    for (label d = 0; d < MeshType::numberOfDirections; d++)
    {
        vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();

        const labelVector pointMeshSize =
            fvMsh_[l].N()
        + MeshType::padding[d]
        + unitXYZ;

        points->SetNumberOfPoints(cmptProduct(pointMeshSize));

        label c = 0;

        for (int k = 0; k < pointMeshSize.z(); k++)
        for (int j = 0; j < pointMeshSize.y(); j++)
        for (int i = 0; i < pointMeshSize.x(); i++)
        {
            // Trim ghost points

            const vector ijk
            (
                i + MeshType::shift[d].x(),
                j + MeshType::shift[d].y(),
                k + MeshType::shift[d].z()
            );

            vector ijkTrim(ijk);

            for (int dim = 0; dim < 3; dim++)
            {
                if (MeshType::shift[d][dim] != 0.0)
                {
                    ijkTrim[dim] =
                        Foam::max
                        (
                            Foam::min
                            (
                                ijkTrim[dim],
                                pointMeshSize[dim]-2.0
                            ),
                            0.0
                        );
                }
            }

            if (ijk != ijkTrim)
            {
                points->SetPoint(c++, p(ijkTrim).v_);
            }
            else
            {
                points->SetPoint(c++, p.interp(i,j,k,MeshType::shift[d]).v_);
            }
        }

        // Create structured grid

        vtkSmartPointer<vtkStructuredGrid> grid =
            vtkSmartPointer<vtkStructuredGrid>::New();

        grid->SetDimensions(pointMeshSize.v_);
        grid->SetPoints(points);

        // Add fields

        writeFields<label,MeshType>(grid.GetPointer(), l, d);
        writeFields<scalar,MeshType>(grid.GetPointer(), l, d);
        writeFields<hexScalar,MeshType>(grid.GetPointer(), l, d);
        writeFields<vector,MeshType>(grid.GetPointer(), l, d);
        writeFields<hexVector,MeshType>(grid.GetPointer(), l, d);
        writeFields<tensor,MeshType>(grid.GetPointer(), l, d);
        writeFields<diagTensor,MeshType>(grid.GetPointer(), l, d);
        writeFields<sphericalTensor,MeshType>(grid.GetPointer(), l, d);
        writeFields<symmTensor,MeshType>(grid.GetPointer(), l, d);

        // Write VTK data

        writeVTKData<MeshType>(timeName, grid.GetPointer(), l, d);
    }
}

template<class MeshType>
void IO::readData(const word timeName, const label l)
{
    for (label d = 0; d < MeshType::numberOfDirections; d++)
    {
        vtkSmartPointer<vtkXMLStructuredGridReader> reader =
            vtkSmartPointer<vtkXMLStructuredGridReader>::New();

        const word vtkFile = vtkFileName<MeshType>(l, d, Pstream::myProcNo());

        const fileName vtkFilePath
        (
            fvMsh_.time().path()/timeName/vtkFile
        );

        reader->SetFileName(vtkFilePath.c_str());
        reader->Update();

        vtkSmartPointer<vtkStructuredGrid> grid = reader->GetOutput();

        readFields<label,MeshType>(grid.GetPointer(),l,d);
        readFields<scalar,MeshType>(grid.GetPointer(),l,d);
        readFields<hexScalar,MeshType>(grid.GetPointer(),l,d);
        readFields<vector,MeshType>(grid.GetPointer(),l,d);
        readFields<hexVector,MeshType>(grid.GetPointer(),l,d);
        readFields<tensor,MeshType>(grid.GetPointer(),l,d);
        readFields<diagTensor,MeshType>(grid.GetPointer(),l,d);
        readFields<sphericalTensor,MeshType>(grid.GetPointer(),l,d);
        readFields<symmTensor,MeshType>(grid.GetPointer(),l,d);
    }
}

template<class Type, class MeshType>
void IO::writeFields
(
    vtkStructuredGrid * grid,
    const label l,
    const label d
) const
{
    typedef meshField<Type,MeshType> FieldType;

    const wordList names(fvMsh_.time().names<FieldType>());

    forAll(names, i)
    {
        const FieldType& field =
            fvMsh_.time().lookupObject<FieldType>(names[i]);

        if (field.writeOpt() == IOobject::AUTO_WRITE)
        {
            writeField(grid, field, l, d);
        }
    }
}

template<class Type, class MeshType>
void IO::readFields
(
    vtkStructuredGrid * grid,
    const label l,
    const label d
) const
{
    typedef meshField<Type, MeshType> FieldType;

    const wordList names(fvMsh_.time().names<FieldType>());

    vtkSmartPointer<vtkCellData> data = grid->GetCellData();

    forAll(names, i)
    {
        FieldType& field =
            fvMsh_.time().lookupObjectRef<FieldType>(names[i]);

        if
        (
            field.readOpt() == IOobject::MUST_READ
         || (
                field.readOpt() == IOobject::READ_IF_PRESENT
             && data->HasArray(names[i].c_str())
            )
        )
        {
            readField(grid, field, l, d);
        }
    }
}

template<class MeshType>
void IO::writeVTKData
(
    const word timeName,
    vtkStructuredGrid * grid,
    const label l,
    const label d
) const
{
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

    const word vtkFile = vtkFileName<MeshType>(l, d, Pstream::myProcNo());

    const fileName vtkFilePath
    (
        fvMsh_.time().path()/timeName/vtkFile
    );

    writer->SetFileName(vtkFilePath.c_str());
    writer->SetInputData(grid);

    if (fvMsh_.time().writeFormat() == IOstream::ASCII)
    {
        writer->SetDataModeToAscii();
    }
    else
    {
        writer->SetDataModeToBinary();
    }

    writer->Write();
}

void IO::readTimeData()
{
    timeData_.clear();

    if (isFile(fvMsh_.time().path()/timeDataFileName))
    {
        IFstream timeDataFile(fvMsh_.time().path()/timeDataFileName);

        timeDataFile >> timeData_;
    }
}

void IO::addToTimeData(const word timeName)
{
    List<Tuple2<scalar,label>> timeData2;

    forAll(timeData_, i)
    {
        const word timeNamei(fvMsh_.time().timeName(timeData_[i].first()));

        if (timeNamei != timeName)
        {
            timeData2.append(timeData_[i]);
        }
    }

    if (timeData2.size() != timeData_.size())
    {
        timeData_ = timeData2;
    }

    timeData_.append
    (
        Tuple2<scalar,label>
        (
            std::stod(timeName),
            Pstream::parRun()
          ? Pstream::nProcs()
          : 1
        )
    );

    sortTimeData();

    // Write to time data file

    OFstream timeDataFile(fvMsh_.time().path()/timeDataFileName);

    timeDataFile << timeData_;
}

void IO::sortTimeData()
{
    SortableList<scalar> times(timeData_.size());

    forAll(times, i)
    {
        times[i] = timeData_[i].first();
    }

    times.sort();

    List<Tuple2<scalar,label>> timeData2(timeData_.size());

    forAll(times, i)
    {
        timeData2[i] = timeData_[times.indices()[i]];
    }

    timeData_ = timeData2;
}

template<class MeshType>
void IO::writeTimeData(const label l)
{
    for (label d = 0; d < MeshType::numberOfDirections; d++)
    {
        const word pvdFile = pvdFileName<MeshType>(l,d);

        OFstream PVDFile(fvMsh_.time().path()/pvdFile);

        PVDFile << "<?xml version=\"1.0\"?>" << endl
                << "<VTKFile type=\"Collection\" version=\"0.1\">" << endl
                << "  <Collection>" << endl;

        forAll(timeData_, i)
        {
            const word timeName(fvMsh_.time().timeName(timeData_[i].first()));

            for (label proc = 0; proc < timeData_[i].second(); proc++)
            {
                PVDFile << "    <DataSet timestep=\"" << timeName << "\" "
                        << "group=\"" << i << "\" part=\"" << proc << "\" "
                        << "file=\"" << timeName << "/"
                        << vtkFileName<MeshType>(l,d,proc) << "\"/>"
                        << endl;
            }
        }

        PVDFile << "  </Collection>" << endl
                << "</VTKFile>" << endl;
    }
}

template<class MeshType>
word IO::vtkFileName
(
    const label l,
    const label d,
    const label proc
) const
{
    return
        dataFileNamePrefix
      + "_"
      + Foam::name(l)
      + "_"
      + MeshType::typeName
      + (
            MeshType::numberOfDirections > 1
          ? word("_" + word(MeshType::directionNames[d]))
          : word("")
        )
      + word("_" + Foam::name(proc) + ".")
      + XMLStructuredGridExtension;
}

template<class MeshType>
word IO::pvdFileName(const label l, const label d) const
{
    return
        dataFileNamePrefix
      + "_"
      + Foam::name(l)
      + "_"
      + MeshType::typeName
      + (
            MeshType::numberOfDirections > 1
          ? word("_" + word(MeshType::directionNames[d]))
          : word("")
        )
      + ".pvd";
}

IO::IO(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh),
    timeData_()
{
    if (Pstream::master())
    {
        readTimeData();
    }
}

IO::~IO()
{}

template<class MeshType>
void IO::write(const label l)
{
    if (fvMsh_.time().writeTime())
    {
        writeNow<MeshType>(l);
    }
}

template<class MeshType>
void IO::writeNow(const label l)
{
    const Time& time = fvMsh_.time();

    const word timeName(time.timeName());

    if (Pstream::master())
    {
        if (!isDir(time.path()/timeName))
        {
            mkDir(time.path()/timeName);
        }
    }

    // Wait

    if (Pstream::parRun())
        returnReduce(1.0,sumOp<scalar>());

    // Write data

    writeData<MeshType>(timeName,l);

    if (Pstream::master())
    {
        addToTimeData(timeName);
        writeTimeData<MeshType>(l);
    }
}

template<class MeshType>
void IO::read(const label l)
{
    const Time& time = fvMsh_.time();

    const word timeName(time.timeName());

    if (Pstream::master())
    {
        if (!isDir(time.path()/timeName) || timeName == "0")
        {
            FatalErrorInFunction
                << "Cannot read from time " << timeName << endl
                << exit(FatalError);
        }
    }

    // Read data

    readData<MeshType>(timeName,l);
}

// Instantiate

template void IO::read<colocated>(const label);
template void IO::read<staggered>(const label);

template void IO::write<colocated>(const label);
template void IO::write<staggered>(const label);

}

}

}
