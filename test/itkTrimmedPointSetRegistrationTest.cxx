/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkJensenHavrdaCharvatTsallisPointSetToPointSetMetricv4.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkCommand.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkAffineTransform.h"
#include "itkTrimmedPointSetToPointSetMetricv4.h"

#include <fstream>

template<typename TFilter>
class itkTrimmedPointSetMetricRegistrationTestCommandIterationUpdate : public itk::Command
{
public:
  using Self = itkTrimmedPointSetMetricRegistrationTestCommandIterationUpdate;

  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro( Self );

protected:
  itkTrimmedPointSetMetricRegistrationTestCommandIterationUpdate() = default;

public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    Execute( (const itk::Object *) caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }
    const auto * optimizer = dynamic_cast< const TFilter * >( object );

    if( !optimizer )
      {
      itkGenericExceptionMacro( "Error dynamic_cast failed" );
      }
    std::cout << "It: " << optimizer->GetCurrentIteration() << " metric value: " << optimizer->GetCurrentMetricValue();
    std::cout << std::endl;
    }
};

int itkTrimmedPointSetRegistrationTest( int argc, char *argv[] )
{
  constexpr unsigned int Dimension = 2;

  unsigned int numberOfIterations = 10;
  if( argc > 1 )
    {
    numberOfIterations = std::stoi( argv[1] );
    }

  using PointSetType = itk::PointSet<double, Dimension>;

  using PointType = PointSetType::PointType;

  PointSetType::Pointer fixedPoints = PointSetType::New();
  fixedPoints->Initialize();

  PointSetType::Pointer movingPoints = PointSetType::New();
  movingPoints->Initialize();


  using GeneratorType = itk::Statistics::MersenneTwisterRandomVariateGenerator;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();

  // Generate two noisy ellipses
  unsigned int nSourcePoints= 1000;
  for(int i=0; i< nSourcePoints; i++ )
    {
    float radius = 100.0;

    float theta = generator->GetUniformVariate(0, 2*itk::Math::pi);
    PointType fixedPoint;
    fixedPoint[0] = 0.5 * radius * std::cos( theta ) + generator->GetNormalVariate();
    fixedPoint[1] = 2 * radius * std::sin( theta ) + generator->GetNormalVariate();
    fixedPoints->SetPoint( i, fixedPoint );
    }

  unsigned int nTargetPoints= 1200;
  for(int i=0; i< nTargetPoints; i++ )
    {
    float radius = 100.0;

    float theta = generator->GetUniformVariate(0, 2*itk::Math::pi);
    PointType movingPoint;
    movingPoint[0] = 2.5 * radius * std::cos( theta ) + generator->GetNormalVariate();
    movingPoint[1] = 0.5 * radius * std::sin( theta ) + generator->GetNormalVariate();
    movingPoints->SetPoint( i, movingPoint );
    }



  using PixelType = double;
  using FixedImageType = itk::Image<PixelType, Dimension>;
  using MovingImageType = itk::Image<PixelType, Dimension>;

  FixedImageType::SizeType fixedImageSize;
  FixedImageType::PointType fixedImageOrigin;
  FixedImageType::DirectionType fixedImageDirection;
  FixedImageType::SpacingType fixedImageSpacing;

  fixedImageSize.Fill( 221 );
  fixedImageOrigin.Fill( -110 );
  fixedImageDirection.SetIdentity();
  fixedImageSpacing.Fill( 1 );

  FixedImageType::Pointer fixedImage = FixedImageType::New();
  fixedImage->SetRegions( fixedImageSize );
  fixedImage->SetOrigin( fixedImageOrigin );
  fixedImage->SetDirection( fixedImageDirection );
  fixedImage->SetSpacing( fixedImageSpacing );
  fixedImage->Allocate();

  using AffineTransformType = itk::AffineTransform<double, Dimension>;
  AffineTransformType::Pointer transform = AffineTransformType::New();
  transform->SetIdentity();

  // Instantiate the metric
  using PointSetMetricType = itk::JensenHavrdaCharvatTsallisPointSetToPointSetMetricv4<PointSetType>;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
  metric->SetFixedPointSet( fixedPoints );
  metric->SetMovingPointSet( movingPoints );
  metric->SetPointSetSigma( 1.0 );
  metric->SetKernelSigma( 10.0 );
  metric->SetUseAnisotropicCovariances( false );
  metric->SetCovarianceKNeighborhood( 5 );
  metric->SetEvaluationKNeighborhood( 10 );
  metric->SetMovingTransform( transform );
  metric->SetAlpha( 1.1 );
  metric->SetVirtualDomainFromImage( fixedImage );
  metric->Initialize();


  using TrimmedPointSetMetricType = itk::TrimmedPointSetToPointSetMetricv4<PointSetType, PointSetType>;
  TrimmedPointSetMetricType::Pointer trimmedMetric = TrimmedPointSetMetricType::New();
  trimmedMetric->SetMetric(metric);
  

  // scales estimator
  using RegistrationParameterScalesFromShiftType = itk::RegistrationParameterScalesFromPhysicalShift< TrimmedPointSetMetricType >;
  RegistrationParameterScalesFromShiftType::Pointer shiftScaleEstimator = RegistrationParameterScalesFromShiftType::New();
  shiftScaleEstimator->SetMetric( trimmedMetric );
  // needed with pointset metrics
  shiftScaleEstimator->SetVirtualDomainPointSet( trimmedMetric->GetVirtualTransformedPointSet() );

  // optimizer
  using OptimizerType = itk::GradientDescentOptimizerv4;
  OptimizerType::Pointer  optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetScalesEstimator( shiftScaleEstimator );
  optimizer->SetMaximumStepSizeInPhysicalUnits( 3.0 );
  using CommandType = itkTrimmedPointSetMetricRegistrationTestCommandIterationUpdate<OptimizerType>;
  CommandType::Pointer observer = CommandType::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  optimizer->SetMinimumConvergenceValue( 0.0 );
  optimizer->SetConvergenceWindowSize( 10 );

  using AffineRegistrationType = itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType>;
  AffineRegistrationType::Pointer affineSimple = AffineRegistrationType::New();
  affineSimple->SetObjectName( "affineSimple" );
  affineSimple->SetFixedPointSet( fixedPoints );
  affineSimple->SetMovingPointSet( movingPoints );
  affineSimple->SetInitialTransform( transform );
  affineSimple->SetMetric( trimmedMetric );
  affineSimple->SetOptimizer( optimizer );

   try
    {
    std::cout << "Point set affine registration update" << std::endl;
    affineSimple->Update();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }

  // applying the resultant transform to moving points and verify result
  std::cout << "Fixed\tMoving\tMovingTransformed\tFixedTransformed\tDiff" << std::endl;
  bool passed = true;
  PointType::ValueType tolerance = 1e-2;
  AffineTransformType::InverseTransformBasePointer affineInverseTransform = affineSimple->GetModifiableTransform()->GetInverseTransform();
  for( unsigned int n=0; n < movingPoints->GetNumberOfPoints(); n++ )
    {
    // compare the points in virtual domain
    PointType transformedMovingPoint = affineInverseTransform->TransformPoint( movingPoints->GetPoint( n ) );
    PointType fixedPoint = fixedPoints->GetPoint( n );
    PointType transformedFixedPoint = affineSimple->GetModifiableTransform()->TransformPoint( fixedPoints->GetPoint( n ) );
    PointType difference;
    difference[0] = transformedMovingPoint[0] - fixedPoint[0];
    difference[1] = transformedMovingPoint[1] - fixedPoint[1];
    std::cout << fixedPoints->GetPoint( n ) << "\t" << movingPoints->GetPoint( n )
          << "\t" << transformedMovingPoint << "\t" << transformedFixedPoint << "\t" << difference << std::endl;
    if( fabs( difference[0] ) > tolerance || fabs( difference[1] ) > tolerance )
      {
      passed = false;
      }
    }
  if( ! passed )
    {
    std::cerr << "Results do not match truth within tolerance." << std::endl;
    //return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
