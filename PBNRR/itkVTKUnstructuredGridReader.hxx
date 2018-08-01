/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVTKUnstructuredGridReader.hxx,v $
  Language:  C++
  Date:      $Date: 2011-10-05 18:01:00 $
  Version:   $Revision: 1.19 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVTKUnstructuredGridReader_hxx
#define __itkVTKUnstructuredGridReader_hxx

#include "itkVTKUnstructuredGridReader.h"
#include <fstream>
#include <stdio.h>
#include <string.h>

namespace itk
{

//
// Constructor
//
template<class TOutputMesh>
VTKUnstructuredGridReader<TOutputMesh>
::VTKUnstructuredGridReader()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer());
}

template<class TOutputMesh>
void
VTKUnstructuredGridReader<TOutputMesh>
::GenerateData()
 {
    typename OutputMeshType::Pointer outputMesh = this->GetOutput();

    outputMesh->SetCellsAllocationMethod(
    OutputMeshType::CellsAllocatedDynamicallyCellByCell );

    if( m_FileName == "" )
    {
        itkExceptionMacro(<< "No input FileName");
    }

    //
    // Read input file
    //
    std::ifstream inputFile( m_FileName.c_str() );

    if( !inputFile.is_open() )
    {
        itkExceptionMacro(<< "Unable to open file\n"
            << "inputFilename= " << m_FileName );
    }

    inputFile.imbue(std::locale::classic());
    std::string line;

    // The first line must be "# vtk DataFile Version x.x" where x.x can
    // vary
    std::getline( inputFile, m_Version, '\n' );
    if (inputFile.fail())
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nUnexpected end-of-file trying to read first line.");
    }
    if (m_Version.find("# vtk DataFile Version ") == std::string::npos)
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nOnly vtk legacy format files can be read."
            << "\nThis file does not start with the line: # vtk DataFile Version x.x where x.x is the version.");
    }

    // Next is a one line description
    std::getline( inputFile, m_Header, '\n' );
    if (inputFile.eof())
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nUnexpected end-of-file trying to read header.");
    }

    // Next is the file format
    std::getline( inputFile, line, '\n' );
    if (inputFile.eof())
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nUnexpected end-of-file trying to file format.");
    }
    if (line.find("ASCII") == std::string::npos)
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nFile format is " << line
            << " but only ASCII files can be read.");
    }

    bool foundPoints = false;
    while( !inputFile.eof() )
    {
        std::getline( inputFile, line, '\n' );

        if( line.find("POINTS") != std::string::npos )
        {
            foundPoints = true;
            break;
        }
    }

    if (!foundPoints)
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nUnexpected end-of-file before finding POINTS.");
    }
    itkDebugMacro("POINTS line" << line );

    std::string pointLine( line, strlen("POINTS "), line.length() );
    itkDebugMacro("pointLine " << pointLine );

    int numberOfPoints = -1;

    if( sscanf(pointLine.c_str(),"%d",&numberOfPoints) != 1 )
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nFailed to read numberOfPoints.\n"
            << "       pointLine= " << pointLine );
    }

    itkDebugMacro("numberOfPoints= " << numberOfPoints );

    if( numberOfPoints < 1 )
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "numberOfPoints < 1"
            << "       numberOfPoints line= " << numberOfPoints );
    }

    outputMesh->GetPoints()->Reserve( numberOfPoints );

    //
    // Load the point coordinates into the itk::Mesh
    //

    PointType point;

    for( int i=0; i < numberOfPoints; i++ )
    {
        inputFile >> point;
        if (inputFile.eof())
        {
            itkExceptionMacro(<< "Error while reading file: " << m_FileName
                << "\nUnexpected end-of-file trying to read points.");
        }
        if (inputFile.fail())
        {
            itkExceptionMacro(<< "Error reading file: " << m_FileName
                << "\nInput could not be interpreted as a point.");
        }
        outputMesh->SetPoint( i, point );
    }

    // Continue searching for the CELLS line
    bool foundCells = false;
    while( !inputFile.eof() )
    {
        std::getline( inputFile, line, '\n' );
        if (line.find("CELLS") != std::string::npos )
        {
            foundCells = true;
            break;
        }
    }

    if (!foundCells)
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nUnexpected end-of-file before finding CELLS.");
    }

    itkDebugMacro( "CELLS line" << line );

    std::string cellsLine( line, strlen("CELLS "), line.length() );
    itkDebugMacro( "cellsLine " << cellsLine);

    //
    // Read the number of cells
    //

    CellIdentifier numberOfCells = 0;
    CellIdentifier numberOfIndices = 0;

    if( sscanf( cellsLine.c_str(), "%ld %ld", &numberOfCells,
        &numberOfIndices ) != 2 )
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nFailed to read numberOfCells from subline2"
            << "\ncellsLine = " << cellsLine );
    }

    itkDebugMacro("numberOfCells " << numberOfCells );
    itkDebugMacro("numberOfIndices " << numberOfIndices );

    if( numberOfCells < 1 )
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nnumberOfCells < 1\nnumberOfCells= "
            << numberOfCells );
    }

    if( numberOfIndices < numberOfCells )
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nnumberOfIndices < numberOfCells\n"
            << "numberOfIndices= " << numberOfIndices << "\n"
            << "numberOfCells= " << numberOfCells );
    }

    //
    // Load the cells into the itk::Mesh
    //

    PointIdentifier numberOfCellPoints;
    long ids[4];

    for(CellIdentifier ii=0; ii < numberOfCells; ii++)
    {
        std::getline( inputFile, line, '\n' );
        if( inputFile.eof() )
        {
            itkExceptionMacro(<< "Error reading file: " << m_FileName
                << "\nFailed to read " << numberOfCells
                << " cells before the end of file."
                << " Only read " << ii+1);
        }

        if( line.find("DATA") != std::string::npos )
        {
            itkExceptionMacro(<< "Error reading file: " << m_FileName
                << "\nRead keyword DATA");
        }

        int got = sscanf( line.c_str(), "%ld %ld %ld %ld %ld", &numberOfCellPoints, &ids[0], &ids[1], &ids[2], &ids[3]);
        if(got != 2 && got != 5)
        {
        	itkExceptionMacro(<< "Error reading file: " << m_FileName
        			<< "\nError parsing VTK_TETRA CELLS or VTK_VERTEX CELLS. Expected 5 or 2 items, respectively, but got "
        			<< got << std::endl
        			<< "Line is: " << line);
        }

        if( numberOfCellPoints != 4 && numberOfCellPoints != 1)
        {
            itkExceptionMacro(<< "Error reading file: " << m_FileName
                << "\nnumberOfCellPoints != 4 and numberOfCellPoints != 1\n"
                << "numberOfCellPoints= " << numberOfCellPoints
                << ". VTKUnstructuredGridReader can only read VTK_TETRA or VTK_VERTEX");
        }

        CellAutoPointer cell;

        // VTK_VERTEX
        if(got == 2)
        {
            if( static_cast<long>(ids[0]) < 0 )
                itkExceptionMacro(<< "Error reading file: " << m_FileName << "point ids must be >= 0.\n" "ids=" << ids[0]);

            if( static_cast<long>(ids[0]) >= numberOfPoints)
                itkExceptionMacro(<< "Error reading file: " << m_FileName << "Point ids must be < number of points: " << numberOfPoints
                     << "\nids= " << ids[0]);

            VertexCellType* vertexCell = new VertexCellType;
            vertexCell->SetPointId( 0, ids[0] );
            cell.TakeOwnership( vertexCell );
        }
        // VTK_TETRA
        else
        {
            if( static_cast<long>(ids[0]) < 0 ||
                static_cast<long>(ids[1]) < 0 ||
                static_cast<long>(ids[2]) < 0 ||
                static_cast<long>(ids[3]) < 0 )
            {
                itkExceptionMacro(<< "Error reading file: " << m_FileName
                    << "point ids must be >= 0.\n"
                    "ids=" << ids[0] << " " << ids[1] << " " << ids[2] << " " << ids[3]);
            }

            if( static_cast<long>(ids[0]) >= numberOfPoints ||
                static_cast<long>(ids[1]) >= numberOfPoints ||
                static_cast<long>(ids[2]) >= numberOfPoints ||
                static_cast<long>(ids[3]) >= numberOfPoints )
            {
                itkExceptionMacro(<< "Error reading file: " << m_FileName
                    << "Point ids must be < number of points: "
                    << numberOfPoints
                    << "\nids= " << ids[0] << " " << ids[1] << " " << ids[2] << " " << ids[3]);
            }
            TetrahedronCellType * tetrahedronCell = new TetrahedronCellType;
            for( PointIdentifier k = 0; k < numberOfCellPoints; k++ )
                tetrahedronCell->SetPointId( k, ids[k] );

            cell.TakeOwnership( tetrahedronCell );
        }

        // Set the cell
        outputMesh->SetCell( ii, cell );
    }


    // Continue searching for the CELL_TYPES line
    bool bfoundCellTypes = false;
    while( !inputFile.eof() )
    {
        std::getline( inputFile, line, '\n' );
        if (line.find("CELL_TYPES") != std::string::npos )
        {
            bfoundCellTypes = true;
            break;
        }
    }

    if (!bfoundCellTypes)
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nUnexpected end-of-file before finding CELL_TYPES.");
    }

    itkDebugMacro( "CELL_TYPES line" << line );

    std::string cellsTypesLine( line, strlen("CELL_TYPES "), line.length() );

    itkDebugMacro( "cellsTypesLine " << cellsTypesLine);


    unsigned int numberOfCellTypes = 0;
    if( sscanf( cellsTypesLine.c_str(), "%d", &numberOfCellTypes) != 1 )
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nFailed to read numberOfCellTypes from subline2"
            << "\ncellsTypesLine = " << cellsTypesLine);
    }


    if( numberOfCellTypes != numberOfCells)
    {
        itkExceptionMacro(<< "Error reading file: " << m_FileName
            << "\nnumberOfCellTypes < numberOfCells ");
    }


    for(CellIdentifier i=0; i<numberOfCellTypes; i++)
    {
        std::getline( inputFile, line, '\n' );
        if( inputFile.eof() )
        {
            itkExceptionMacro(<< "Error reading file: " << m_FileName
                << "\nFailed to read " << numberOfCellTypes
                << " cells before the end of file."
                << " Only read " << i+1);
        }

        int typeId=0;
        int got;
        if( (got = sscanf( line.c_str(), "%d ", &typeId )) != 1 )
        {
            itkExceptionMacro(<< "Error reading file: " << m_FileName
                << "\nError parsing Cells types. Expected 1 item but got "
                << got << std::endl
                << "Line is: " << line);
        }
        if( typeId != 10 && typeId != 1)
        {
            itkExceptionMacro(<< "Error reading file: " << m_FileName
                << "\nCell type != 10 and Cell type != 1\n"
                << "Cell type = " << typeId
                << ". VTKUnstructuredGridReader can only read VTK_TETRA or VTK_VERTEX");
        }
    }

    bool foundPointData = false;
    bool foundCellData = false;

    while( !inputFile.eof() )
    {
        std::getline( inputFile, line, '\n' );

        if( line.find("POINT_DATA") != std::string::npos )
        {
            foundPointData = true;
            break;
        }
        if( line.find("CELL_DATA") != std::string::npos )
        {
            foundCellData = true;
            break;
        }
    }

    if( foundPointData )
    {
        typedef typename OutputMeshType::PointDataContainer MyPointDataContainer;

        outputMesh->SetPointData( MyPointDataContainer::New() );
        outputMesh->GetPointData()->Reserve( numberOfPoints );

        itkDebugMacro("POINT_DATA line" << line );

        // Skip two lines
        for(int j=0; j<2; j++)
        {
            if (!inputFile.eof())
            {
                std::getline( inputFile, line, '\n' );
            }
            else
            {
                //itkExceptionMacro(<< "Error reading file: " << m_FileName
                //<< "\nUnexpected end-of-file while trying to read POINT_DATA.");
                inputFile.close();
                return;
            }
        }


        double pointData;

        for( int pid=0; pid < numberOfPoints; pid++ )
        {
            if (inputFile.eof())
            {
                break;
                //itkExceptionMacro(<< "Error reading file: " << m_FileName
                //<< "\nUnexpected end-of-file while trying to read POINT_DATA."
                //<< "Failed while trying to reading point data for id: " << pid);
            }
            inputFile >> pointData;
            outputMesh->SetPointData( pid, pointData );
        }
    }

    if( foundCellData )
    {
        typedef typename OutputMeshType::CellDataContainer MyCellDataContainer;

        outputMesh->SetCellData( MyCellDataContainer::New() );
        outputMesh->GetCellData()->Reserve( numberOfCells);

        itkDebugMacro("CELL_DATA line" << line );

        // Skip two lines
        for(int j=0; j<2; j++)
        {
            if (!inputFile.eof())
            {
                std::getline( inputFile, line, '\n' );
            }
            else
            {
                //itkExceptionMacro(<< "Error reading file: " << m_FileName
                //<< "\nUnexpected end-of-file while trying to read POINT_DATA.");
                inputFile.close();
                return;
            }
        }

        double cellData;

        for( int pid=0; pid < (int)numberOfCells; pid++ )
        {
            if (inputFile.eof()){
                break;
                //itkExceptionMacro(<< "Error reading file: " << m_FileName
                //<< "\nUnexpected end-of-file while trying to read POINT_DATA."
                //<< "Failed while trying to reading point data for id: " << pid);
            }
            inputFile >> cellData;
            outputMesh->SetCellData( pid, cellData );
        }
    }

    inputFile.close();
 }

template<class TOutputMesh>
void
VTKUnstructuredGridReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " << m_FileName << std::endl;
  os << indent << "Version: " << m_Version << std::endl;
  os << indent << "Header: " << m_Header << std::endl;
}

} //end of namespace itk


#endif
