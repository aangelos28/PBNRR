#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "PBNRRCLP.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMesh.h"
#include "itkVTKUnstructuredGridReader.h"
#include "itkPhysicsBasedNonRigidRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"

#include "vtkUnstructuredGrid.h"
#include "vtkIdList.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkNew.h"


// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

void PrintLinearTetrahedraVTK(itk::VectorContainer < unsigned long ,itk::fem::Element::Node::Pointer> *NodeContainer,
		itk::VectorContainer<unsigned long  ,itk::fem::Element::Pointer> *ElementContainer,std::string OutputFileName)
{
	if(NodeContainer == NULL || ElementContainer == NULL || OutputFileName.empty())
		return;

	vtkNew<vtkUnstructuredGrid> pUG;
	vtkPoints* pPoints = vtkPoints::New();

	// Set nodes coordinates
	int numNodes = NodeContainer->Size();

	itk::fem::Element::Node::Pointer pNode = NULL;
	itk::fem::Element::VectorType vec;

    for( int i=0; i < numNodes; i++ )
    {
    	pNode = NodeContainer->at(i);
    	vec  =  pNode->GetCoordinates();
    	pPoints->InsertPoint(i,vec[0],vec[1],vec[2]);
    }

    // Set the points to Unstructured Grid
	pUG->SetPoints(pPoints);
	pPoints->Delete();

	// Set the id's of cells vertices
	vtkIdList* ptIds = vtkIdList::New();

	itk::fem::Element::Pointer pTetr = NULL;
	long int ids[4];

	int numTets = ElementContainer->Size();

	for( int i=0; i < numTets; i++ )
	{
		pTetr = ElementContainer->at(i);

		ids[0] = pTetr->GetNode(0)->GetGlobalNumber();
		ids[1] = pTetr->GetNode(1)->GetGlobalNumber();
		ids[2] = pTetr->GetNode(2)->GetGlobalNumber();
		ids[3] = pTetr->GetNode(3)->GetGlobalNumber();

		ptIds->InsertId(0,ids[0]);
		ptIds->InsertId(1,ids[1]);
		ptIds->InsertId(2,ids[2]);
		ptIds->InsertId(3,ids[3]);

		// Insert the cell to Unstructured Grid (10 = VTK_TETRA)
		pUG->InsertNextCell(10,ptIds);
	}

	ptIds->Delete();

	vtkUnstructuredGridWriter* pWriter = vtkUnstructuredGridWriter::New();
	pWriter->SetFileTypeToASCII();
	pWriter->SetFileName(OutputFileName.c_str());
	pWriter->SetInputData(pUG.GetPointer());
	pWriter->Write();
	pWriter->Update();
	pWriter->Delete();
}

template <typename TFixedImage, typename TMovingImage, typename TMaskImage, typename TMesh, typename TDeformationField>
typename TDeformationField::Pointer CreateInverseDeformationField(typename itk::fem::PhysicsBasedNonRigidRegistrationMethod<TFixedImage, TMovingImage, TMaskImage, TMesh, TDeformationField>::Pointer pPBNRRFilter)
{
	if( pPBNRRFilter.IsNull() )
		return NULL;

	typename TDeformationField::Pointer pDirectField = pPBNRRFilter->GetOutput();
	if(pDirectField.IsNull())
		return NULL;

	const unsigned int ImageDimension = 3;
	itk::fem::Element::VectorType vMin,vMax,VoxelDispl,nodeDispl,vGlobal,vLocal,shapeF,dVec;
	itk::fem::Element::Pointer elem;
	typename TDeformationField::RegionType region;
	typename TDeformationField::SizeType   regionSize;
	typename TDeformationField::PixelType  PixelValue;
	std::vector<itk::fem::Element::VectorType> nodeDisplVector;
	unsigned int numDofs = 0;
	itk::fem::Element::DegreeOfFreedomIDType id;
	typename TDeformationField::PointType global_pt;
	typename TDeformationField::PointType point;

	// Allocate space for the warped image
	typename TDeformationField::Pointer pInverseField = TDeformationField::New();
	region.SetIndex(pDirectField->GetLargestPossibleRegion().GetIndex());
	region.SetSize(pDirectField->GetLargestPossibleRegion().GetSize());
	pInverseField->SetRegions(region);
	pInverseField->SetOrigin(pDirectField->GetOrigin());
	pInverseField->SetSpacing(pDirectField->GetSpacing());
	pInverseField->SetDirection(pDirectField->GetDirection());
	pInverseField->Allocate();
	PixelValue.Fill(0);
	pInverseField->FillBuffer(PixelValue);

	typename itk::fem::FEMObject<ImageDimension>::ElementContainerType::Iterator it;
	typename itk::fem::PhysicsBasedNonRigidRegistrationMethod<TFixedImage, TMovingImage, TMaskImage, TMesh, TDeformationField>::FEMFilterType::FEMSolverType::FEMObjectType* pFEMObject =
			pPBNRRFilter->GetFEMFilter()->GetFEMSolver()->GetOutput();

	if(!pFEMObject)
		return NULL;

	// Loop through all the cells of mesh
	for(it = pFEMObject->GetElementContainer()->Begin(); it != pFEMObject->GetElementContainer()->End(); ++it)
	{
		elem = it.Value();

		// Initialize the bounding box
		vMin = elem->GetNodeCoordinates(0);
		vMax = vMin;
		const unsigned int NumberOfDimensions = elem->GetNumberOfSpatialDimensions();

		// Find cell's bounding box
		for( unsigned int i = 1; i < elem->GetNumberOfNodes(); i++ )
		{
			const itk::fem::Element::VectorType& v = elem->GetNodeCoordinates(i);
			for( unsigned int j = 0; j < NumberOfDimensions; j++ )
			{
				if( v[j] < vMin[j] )
					vMin[j] = v[j];
				if( v[j] > vMax[j] )
					vMax[j] = v[j];
			}
		}

		numDofs = elem->GetNumberOfDegreesOfFreedomPerNode();
		nodeDispl.set_size(numDofs);
		nodeDisplVector.clear();

		// Take the displacements of each node
		for( unsigned int i = 0; i < elem->GetNumberOfNodes(); i++ )
		{
			nodeDispl.fill(0);
			for(int j=0; j<(int)numDofs; j++)
			{
				id = elem->GetNode(i)->GetDegreeOfFreedom(j);
				nodeDispl[j] = pPBNRRFilter->GetFEMFilter()->GetFEMSolver()->GetSolution(id,0);
			}
			nodeDisplVector.push_back(nodeDispl);
		}
		typename TDeformationField::PointType minPoint,maxPoint;
		for( unsigned int j = 0; j < ImageDimension; j++ )
		{
			minPoint[j] = vMin[j];
			maxPoint[j] = vMax[j];
		}

		// Check if the two corners of the bounding box is inside the moving image.
		typename TDeformationField::IndexType vMinIndex,vMaxIndex,vDirectionIndex;
		if( !pDirectField->TransformPhysicalPointToIndex(minPoint,vMinIndex) )
			continue;
		if( !pDirectField->TransformPhysicalPointToIndex(maxPoint,vMaxIndex) )
			continue;

		// Set the region size of the bounding box
		for( unsigned int i = 0; i < NumberOfDimensions; i++ )
			regionSize[i] = vMaxIndex[i] - vMinIndex[i] + 1;

		for( unsigned int i = 0; i < NumberOfDimensions; i++ )
			vDirectionIndex[i] = vMinIndex[i];

		region.SetSize(regionSize);
		region.SetIndex(vDirectionIndex);

		itk::ImageRegionIterator<TDeformationField> iter(pInverseField,region);
		typename TDeformationField::IndexType index;
		vGlobal.set_size(numDofs);
		dVec.set_size(numDofs);
		dVec.fill(0);

		// Step over all voxels of the region size
		for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
		{
			index = iter.GetIndex();
			global_pt.Fill(0);

			// Get the physical coordinates of the voxel
			pInverseField->TransformIndexToPhysicalPoint(index,global_pt);
			for( unsigned int j = 0; j < ImageDimension; j++ )
				vGlobal[j] = global_pt[j];

			// Check if the point is within the element...
			if(elem->GetLocalFromGlobalCoordinates(vGlobal,vLocal))
			{
				// get the shape functions
				shapeF.fill(0);
				shapeF = elem->ShapeFunctions(vLocal);
				dVec.fill(0);

				// u = u1*L1 + u2*L2 + u3*L3 + u4*L4
				for(int k=0; k<(int)elem->GetNumberOfNodes(); k++)
					dVec = dVec + (shapeF[k] * nodeDisplVector.at(k));

	            for(unsigned int i = 0; i < ImageDimension; i++)
	            	PixelValue[i] = - dVec[i];

	            pInverseField->SetPixel(index, PixelValue);
			}
		}
	}
	return pInverseField;
}


template <typename TFixedImage, typename TMovingImage, typename TMaskImage, typename TMesh, typename TDeformationField>
typename TMovingImage::Pointer CreateWarpedImage(typename itk::fem::PhysicsBasedNonRigidRegistrationMethod<TFixedImage, TMovingImage, TMaskImage, TMesh, TDeformationField>::Pointer pPBNRRFilter)
{
	if( pPBNRRFilter.IsNull() )
		return NULL;

	const TMovingImage* pMovingImage = pPBNRRFilter->GetMovingImage();
	if(!pMovingImage)
		return NULL;

	const unsigned int ImageDimension = 3;
	itk::fem::Element::VectorType vMin,vMax,VoxelDispl,nodeDispl,vGlobal,vLocal,shapeF,dVec,pos;
	itk::fem::Element::Pointer elem;
	typename TMovingImage::RegionType region;
	typename TMovingImage::SizeType   regionSize;
	typename TMovingImage::PixelType  PixelValue;
	std::vector<itk::fem::Element::VectorType> nodeDisplVector;
	unsigned int numDofs = 0;
	itk::fem::Element::DegreeOfFreedomIDType id;
	typename TMovingImage::PointType global_pt;
	typename TMovingImage::PointType point;

	// Allocate space for the warped image
	typename TMovingImage::Pointer pWarpedImage = TMovingImage::New();
	region.SetIndex(pMovingImage->GetLargestPossibleRegion().GetIndex());
	region.SetSize(pMovingImage->GetLargestPossibleRegion().GetSize());
	pWarpedImage->SetRegions(region);
	pWarpedImage->SetOrigin(pMovingImage->GetOrigin());
	pWarpedImage->SetSpacing(pMovingImage->GetSpacing());
	pWarpedImage->SetDirection(pMovingImage->GetDirection());
	pWarpedImage->Allocate();
	pWarpedImage->FillBuffer(0);

	// Set the interpolator as linear
	typedef itk::LinearInterpolateImageFunction<TMovingImage, double> InterpolatorType;
	typename InterpolatorType::Pointer pInterpolator = InterpolatorType::New();
	pInterpolator->SetInputImage(pMovingImage);

	typename itk::fem::FEMObject<ImageDimension>::ElementContainerType::Iterator it;
	typename itk::fem::PhysicsBasedNonRigidRegistrationMethod<TFixedImage, TMovingImage, TMaskImage, TMesh, TDeformationField>::FEMFilterType::FEMSolverType::FEMObjectType* pFEMObject =
			pPBNRRFilter->GetFEMFilter()->GetFEMSolver()->GetOutput();

	if(!pFEMObject)
		return NULL;

	// Loop through all the cells of mesh
	for(it = pFEMObject->GetElementContainer()->Begin(); it != pFEMObject->GetElementContainer()->End(); ++it)
	{
		elem = it.Value();

		// Initialize the bounding box
		vMin = elem->GetNodeCoordinates(0);
		vMax = vMin;
		const unsigned int NumberOfDimensions = elem->GetNumberOfSpatialDimensions();

		// Find cell's bounding box
		for( unsigned int i = 1; i < elem->GetNumberOfNodes(); i++ )
		{
			const itk::fem::Element::VectorType& v = elem->GetNodeCoordinates(i);
			for( unsigned int j = 0; j < NumberOfDimensions; j++ )
			{
				if( v[j] < vMin[j] )
					vMin[j] = v[j];
				if( v[j] > vMax[j] )
					vMax[j] = v[j];
			}
		}

		numDofs = elem->GetNumberOfDegreesOfFreedomPerNode();
		nodeDispl.set_size(numDofs);
		nodeDisplVector.clear();

		// Take the displacements of each node
		for( unsigned int i = 0; i < elem->GetNumberOfNodes(); i++ )
		{
			nodeDispl.fill(0);
			for(int j=0; j<(int)numDofs; j++)
			{
				id = elem->GetNode(i)->GetDegreeOfFreedom(j);
				nodeDispl[j] = pPBNRRFilter->GetFEMFilter()->GetFEMSolver()->GetSolution(id,0);
			}
			nodeDisplVector.push_back(nodeDispl);
		}
		typename TMovingImage::PointType minPoint,maxPoint;
		for( unsigned int j = 0; j < ImageDimension; j++ )
		{
			minPoint[j] = vMin[j];
			maxPoint[j] = vMax[j];
		}

		// Check if the two corners of the bounding box is inside the moving image.
		typename TMovingImage::IndexType vMinIndex,vMaxIndex,vDirectionIndex;
		if( !pMovingImage->TransformPhysicalPointToIndex(minPoint,vMinIndex) )
			continue;
		if( !pMovingImage->TransformPhysicalPointToIndex(maxPoint,vMaxIndex) )
			continue;

		// Set the region size of the bounding box
		for( unsigned int i = 0; i < NumberOfDimensions; i++ )
			regionSize[i] = vMaxIndex[i] - vMinIndex[i] + 1;

		for( unsigned int i = 0; i < NumberOfDimensions; i++ )
			vDirectionIndex[i] = vMinIndex[i];

		region.SetSize(regionSize);
		region.SetIndex(vDirectionIndex);

		itk::ImageRegionIterator<TMovingImage> iter((TMovingImage*)pPBNRRFilter->GetMovingImage(),region);
		typename TMovingImage::IndexType index;
		vGlobal.set_size(numDofs);
		dVec.set_size(numDofs);
		dVec.fill(0);

		// Step over all voxels of the region size
		for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
		{
			index = iter.GetIndex();
			global_pt.Fill(0);

			// Get the physical coordinates of the voxel
			pMovingImage->TransformIndexToPhysicalPoint(index,global_pt);
			for( unsigned int j = 0; j < ImageDimension; j++ )
				vGlobal[j] = global_pt[j];

			// Check if the point is within the element...
			if(elem->GetLocalFromGlobalCoordinates(vGlobal,vLocal))
			{
				// get the shape functions
				shapeF.fill(0);
				shapeF = elem->ShapeFunctions(vLocal);
				dVec.fill(0);

				// u = u1*L1 + u2*L2 + u3*L3 + u4*L4
				for(int k=0; k<(int)elem->GetNumberOfNodes(); k++)
					dVec = dVec + (shapeF[k] * nodeDisplVector.at(k));

				pos = vGlobal - dVec;

				for(int k=0; k<(int)ImageDimension; k++)
					point[k] = pos[k];

				// Point inside the image buffer
				if( pInterpolator->IsInsideBuffer(point) )
				{
					PixelValue = pInterpolator->Evaluate(point);
					pWarpedImage->SetPixel(index, PixelValue);
				}
			}
		}
	}
	return pWarpedImage;
}

template<class TPixel, unsigned int  VImageDimension >
bool WriteImage(typename itk::Image< TPixel, VImageDimension >::Pointer pImage ,std::string fname)
{
	typedef itk::ImageFileWriter< itk::Image< TPixel,VImageDimension > > WriterType;
	typename WriterType::Pointer pWriter = WriterType::New();
	pWriter->SetFileName(fname);
	pWriter->SetInput(pImage);
	try
	{
		pWriter->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Error during writing image: " << e << std::endl;
		return false;
	}
	return true;
}

template <class T1, class T2, class T3>
int DoIt3( int argc, char * argv[], const T1 &, const T2 &, const T3 & )
{
	//
	// Command line processing
	//
	PARSE_ARGS;

	const    unsigned int ImageDimension = 3;
	typedef  T1                                        				FixedPixelType;
	typedef  T2                                        				MovingPixelType;
	typedef  T3                                        				MaskPixelType;
	typedef float 													MeshPixelType;
	typedef itk::Vector< MeshPixelType, ImageDimension > 			DeformationFieldPixelType;

	typedef itk::Image<FixedPixelType, ImageDimension> 				FixedImageType;
	typedef itk::Image<MovingPixelType, ImageDimension> 			MovingImageType;
	typedef itk::Image<MaskPixelType, ImageDimension> 				MaskImageType;
	typedef itk::Image< DeformationFieldPixelType,ImageDimension >  DeformationFieldType;
	typedef itk::Mesh<MeshPixelType,ImageDimension> 				MeshType;

	typedef itk::ImageFileReader<FixedImageType>          			FixedImageReaderType;
	typedef itk::ImageFileReader<MovingImageType>         			MovingImageReaderType;
	typedef itk::ImageFileReader<MaskImageType>           			MaskImageReaderType;


	// Add a time probe
	itk::TimeProbesCollectorBase collector;

	// Read the fixed volume
	typename FixedImageReaderType::Pointer pFixedReader = FixedImageReaderType::New();
	pFixedReader->SetFileName( FixedImageFileName.c_str() );
	try
	{
	    collector.Start( "Read fixed image" );
		pFixedReader->Update();
	    collector.Stop( "Read fixed image" );
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error Reading Fixed image: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	typename FixedImageType::Pointer pFixedImage = pFixedReader->GetOutput();

	// Read the moving volume
	typename MovingImageReaderType::Pointer pMovingReader = MovingImageReaderType::New();
	pMovingReader->SetFileName( MovingImageFileName.c_str() );
	try
	{
	    collector.Start( "Read moving image" );
		pMovingReader->Update();
	    collector.Stop( "Read moving image" );
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error Reading Moving image: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	typename MovingImageType::Pointer pMovingImage = pMovingReader->GetOutput();

	// Read the mask volume
	typename MaskImageReaderType::Pointer pMaskReader = MaskImageReaderType::New();
	pMaskReader->SetFileName( MaskImageFileName.c_str() );
	try
	{
	    collector.Start( "Read mask image" );
		pMaskReader->Update();
	    collector.Stop( "Read mask image" );
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error Reading Mask image: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	typename MaskImageType::Pointer pMaskImage = pMaskReader->GetOutput();

	// Read the Tetrahedral mesh
	typename itk::VTKUnstructuredGridReader<MeshType>::Pointer pMeshReader = itk::VTKUnstructuredGridReader<MeshType>::New();
	pMeshReader->SetFileName(MeshFileName);
	try
	{
		collector.Start( "Read mesh" );
		pMeshReader->Update();
		collector.Stop( "Read mesh" );
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "Error Reading Mesh: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	MeshType::Pointer pMesh = pMeshReader->GetOutput();

	std::cout << "*********************************************" << std::endl;
	std::cout << "Fixed Image : " << FixedImageFileName << std::endl;
	std::cout << "Moving Image : " << MovingImageFileName << std::endl;
	std::cout << "Mask Image : " << MaskImageFileName << std::endl;
	std::cout << "Mesh : " << MeshFileName << std::endl;
	std::cout << "Young's Modulus : " << YoungModulus << std::endl;
	std::cout << "Poisson's Ratio : " << PoissonRatio << std::endl;
    std::cout << "Block Radius X : " << BlockRadius[0] << std::endl;
    std::cout << "Block Radius Y : " << BlockRadius[1] << std::endl;
    std::cout << "Block Radius Z : " << BlockRadius[2] << std::endl;
    std::cout << "Search Radius X : " << WindowRadius[0] << std::endl;
    std::cout << "Search Radius Y : " << WindowRadius[1] << std::endl;
    std::cout << "Search Radius Z : " << WindowRadius[2] << std::endl;
    std::cout << "Selection Fraction : " << SelectionFraction << std::endl;
    std::cout << "Rejection Fraction : " << RejectionFraction << std::endl;
    std::cout << "Number Outlier Rejection Steps : " << OutlierRejectionSteps << std::endl;
    std::cout << "Number Interpolation Steps : " << InterpolationSteps << std::endl;
	std::cout << "Warped Image : " << WarpedMovingImageFileName << std::endl;
	std::cout << "Direct Deformation Field : " << DirectDeformationFieldFileName << std::endl;
	std::cout << "Inverse Deformation Field : " << InverseDeformationFieldFileName << std::endl;
	std::cout << "Warped Mesh : " << WarpedMeshFileName << std::endl;
	std::cout << "*********************************************" << std::endl;
	std::cout << " " << std::endl;

	// Check parameters
	if(YoungModulus < 1e-6)
	{
	    std::cerr << "Invalid value for YoungModulus :  " << YoungModulus
	              << "\nShould be greater than 0." << std::endl;
	    return EXIT_FAILURE;
	}
	if(PoissonRatio < 0.1 || PoissonRatio > 0.49 )
	{
	    std::cerr << "Invalid value for PoissonRatio :  " << PoissonRatio
	              << "\nShould be greater than 0.1 and lower than 0.49" << std::endl;
	    return EXIT_FAILURE;
	}
	if(BlockRadius[0] < 1 || BlockRadius[1] < 1 || BlockRadius[2] < 1)
	{
	    std::cerr << "Invalid value for BlockRadius :  " << BlockRadius[0] << "," << BlockRadius[1] << "," << BlockRadius[2]
	              << "\nShould be at least 1,1,1" << std::endl;
	    return EXIT_FAILURE;
	}
	if(WindowRadius[0] < 1 || WindowRadius[1] < 1 || WindowRadius[2] < 1)
	{
	    std::cerr << "Invalid value for WindowRadius :  " << WindowRadius[0] << "," << WindowRadius[1] << "," << WindowRadius[2]
	              << "\nShould be at least 1,1,1" << std::endl;
	    return EXIT_FAILURE;
	}
	if(SelectionFraction < 0.01 || SelectionFraction > 1.0 )
	{
	    std::cerr << "Invalid value for SelectionFraction :  " << SelectionFraction
	              << "\nShould be greater or equal to 0.01 and lower or equal to 1." << std::endl;
	    return EXIT_FAILURE;
	}
	if(RejectionFraction < 0.01 || RejectionFraction >= 1.0 )
	{
	    std::cerr << "Invalid value for RejectionFraction :  " << RejectionFraction
	              << "\nShould be greater or equal to 0.01 and lower than 1." << std::endl;
	    return EXIT_FAILURE;
	}
	if(OutlierRejectionSteps < 1)
	{
	    std::cerr << "Invalid value for OutlierRejectionSteps :  " << OutlierRejectionSteps
	              << "\nShould be greater than 0." << std::endl;
	    return EXIT_FAILURE;
	}
	if(InterpolationSteps < 0)
	{
	    std::cerr << "Invalid value for InterpolationSteps :  " << InterpolationSteps
	              << "\nShould be non negative." << std::endl;
	    return EXIT_FAILURE;
	}

	// Create the filter
	typedef itk::fem::PhysicsBasedNonRigidRegistrationMethod<FixedImageType, MovingImageType, MaskImageType, MeshType, DeformationFieldType>  PBNRRFilterType;
	typename PBNRRFilterType::Pointer pNRRFilter = PBNRRFilterType::New();

	pNRRFilter->SetFixedImage(pFixedImage);
	pNRRFilter->SetMovingImage(pMovingImage);
	pNRRFilter->SetMaskImage(pMaskImage);
	pNRRFilter->SetMesh(pMesh);
	pNRRFilter->SetSelectFraction(SelectionFraction);
	pNRRFilter->SetRejectFraction(RejectionFraction);
	pNRRFilter->SetYoungModulus(YoungModulus);
	pNRRFilter->SetPoissonRatio(PoissonRatio);

	itk::Size< ImageDimension > blockHalfSize;
	itk::Size< ImageDimension > windowHalfSize;
	for(unsigned int i=0; i<ImageDimension; i++)
	{
		blockHalfSize[i] = BlockRadius[i];
		windowHalfSize[i] = WindowRadius[i];
	}

	pNRRFilter->SetBlockRadius(blockHalfSize);
	pNRRFilter->SetSearchRadius(windowHalfSize);
	pNRRFilter->SetOutlierRejectionSteps(OutlierRejectionSteps);
	pNRRFilter->SetApproximationSteps(InterpolationSteps);

	std::cout << "Filter: " << pNRRFilter << std::endl;
	try
	{
	    collector.Start( "Register" );
		pNRRFilter->Update();
	    collector.Stop( "Register" );
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "Error Update the NRR filter: " << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	// Write the warped image
	if(!WarpedMovingImageFileName.empty())
	{
		std::cout << "Creating the Warped Image..." << std::endl;
		collector.Start( "Create Warped image" );
		typename MovingImageType::Pointer pWarpedImage = CreateWarpedImage<FixedImageType, MovingImageType, MaskImageType, MeshType, DeformationFieldType>(pNRRFilter);
		collector.Stop( "Create Warped image" );

		std::cout << "Write Warped Image..." << std::endl;

		collector.Start( "Write Warped image" );
		if(!WriteImage<MovingPixelType,ImageDimension>(pWarpedImage,WarpedMovingImageFileName))
			return EXIT_FAILURE;
		collector.Stop( "Write Warped image" );

	}

	// Write the direct deformation field
	if(!DirectDeformationFieldFileName.empty())
	{
		std::cout << "Writing Direct Deformation Field..." << std::endl;
		collector.Start( "Write Direct Image Deformation Field " );
		if(!WriteImage<DeformationFieldPixelType,ImageDimension>(pNRRFilter->GetOutput(),DirectDeformationFieldFileName))
			return EXIT_FAILURE;
		collector.Stop( "Write Direct Image Deformation Field " );
	}

	// Write the inverse deformation field
	if(!InverseDeformationFieldFileName.empty())
	{
		std::cout << "Creating the Inverse Deformation Field..." << std::endl;
		collector.Start( "Compute Inverse Image Deformation Field " );
		typename DeformationFieldType::Pointer pInverseField = CreateInverseDeformationField<FixedImageType, MovingImageType, MaskImageType, MeshType, DeformationFieldType>(pNRRFilter);
		collector.Stop( "Compute Inverse Image Deformation Field " );

		std::cout << "Writing the Inverse Deformation Field..." << std::endl;
		collector.Start( "Write Inverse Image Deformation Field " );
		if(!WriteImage<DeformationFieldPixelType,ImageDimension>(pInverseField,InverseDeformationFieldFileName))
			return EXIT_FAILURE;
		collector.Stop( "Write Inverse Image Deformation Field " );
	}

	// Write the warped mesh
	if(!WarpedMeshFileName.empty())
	{
		std::cout << "Writing the warped mesh..." << std::endl;
		collector.Start( "Write warped mesh " );
		PrintLinearTetrahedraVTK(pNRRFilter->GetFEMFilter()->GetFEMSolver()->GetOutput()->GetNodeContainer(),
				pNRRFilter->GetFEMFilter()->GetFEMSolver()->GetOutput()->GetElementContainer(),
				WarpedMeshFileName);
		collector.Stop( "Write warped mesh " );
	}

	// Report the time taken by the registration
	collector.Report();

	return EXIT_SUCCESS;
}

template <class T1, class T2>
int DoIt2( int argc, char * argv[], const T1 & targ1, const T2 & targ2)
{
	PARSE_ARGS;

	itk::ImageIOBase::IOPixelType     pixelType;
	itk::ImageIOBase::IOComponentType componentType;

	try
	{
		itk::GetImageType(MaskImageFileName, pixelType, componentType);

		// This filter handles all types

		switch( componentType )
		{
		case itk::ImageIOBase::CHAR:
		case itk::ImageIOBase::UCHAR:
		case itk::ImageIOBase::SHORT:
			return DoIt3( argc, argv, targ1, targ2, static_cast<short>(0) );
			break;
		case itk::ImageIOBase::USHORT:
		case itk::ImageIOBase::INT:
			return DoIt3( argc, argv, targ1, targ2, static_cast<int>(0) );
			break;
		case itk::ImageIOBase::UINT:
		case itk::ImageIOBase::ULONG:
			return DoIt3( argc, argv, targ1, targ2, static_cast<unsigned long>(0) );
			break;
		case itk::ImageIOBase::LONG:
			return DoIt3( argc, argv, targ1, targ2, static_cast<long>(0) );
			break;
		case itk::ImageIOBase::FLOAT:
			return DoIt3( argc, argv, targ1, targ2, static_cast<float>(0) );
			break;
		case itk::ImageIOBase::DOUBLE:
			return DoIt3( argc, argv, targ1, targ2, static_cast<float>(0) );
			break;
		case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
		default:
			std::cout << "unknown component type" << std::endl;
			break;
		}
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << argv[0] << ": exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_FAILURE;
}

template <class T>
int DoIt( int argc, char * argv[], const T& targ)
{

  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(MovingImageFileName, pixelType, componentType);

    // This filter handles all types

    switch( componentType )
      {
      case itk::ImageIOBase::CHAR:
      case itk::ImageIOBase::UCHAR:
      case itk::ImageIOBase::SHORT:
        return DoIt2( argc, argv, targ, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::USHORT:
      case itk::ImageIOBase::INT:
        return DoIt2( argc, argv, targ, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::UINT:
      case itk::ImageIOBase::ULONG:
        return DoIt2( argc, argv, targ, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt2( argc, argv, targ, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt2( argc, argv, targ, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt2( argc, argv, targ, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_FAILURE;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
	PARSE_ARGS;

	itk::ImageIOBase::IOPixelType     pixelType;
	itk::ImageIOBase::IOComponentType componentType;

	try
	{
		itk::GetImageType(FixedImageFileName, pixelType, componentType);

		// This filter handles all types

		switch( componentType )
		{
		case itk::ImageIOBase::CHAR:
		case itk::ImageIOBase::UCHAR:
		case itk::ImageIOBase::SHORT:
			return DoIt( argc, argv, static_cast<short>(0) );
			break;
		case itk::ImageIOBase::USHORT:
		case itk::ImageIOBase::INT:
			return DoIt( argc, argv, static_cast<int>(0) );
			break;
		case itk::ImageIOBase::UINT:
		case itk::ImageIOBase::ULONG:
			return DoIt( argc, argv, static_cast<unsigned long>(0) );
			break;
		case itk::ImageIOBase::LONG:
			return DoIt( argc, argv, static_cast<long>(0) );
			break;
		case itk::ImageIOBase::FLOAT:
			return DoIt( argc, argv, static_cast<float>(0) );
			break;
		case itk::ImageIOBase::DOUBLE:
			return DoIt( argc, argv, static_cast<float>(0) );
			break;
		case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
		default:
			std::cout << "unknown component type" << std::endl;
			break;
		}
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << argv[0] << ": exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
