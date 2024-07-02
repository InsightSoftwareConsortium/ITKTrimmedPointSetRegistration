/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkTrimmedEuclideanDistancePointSetToPointSetMetricv4_h
#define itkTrimmedEuclideanDistancePointSetToPointSetMetricv4_h

#include "itkEuclideanDistancePointSetToPointSetMetricv4.h"

#include "itkTimeProbe.h"

namespace itk
{
/** \class TrimmedEuclideanDistancePointSetToPointSetMetricv4
 * \brief Computes similarity between two point sets.
 *
 * \ingroup TrimmedPointSetRegistration
 */

template<typename TFixedPointSet,  typename TMovingPointSet = TFixedPointSet,
  class TInternalComputationValueType = double>
class ITK_TEMPLATE_EXPORT TrimmedEuclideanDistancePointSetToPointSetMetricv4
: public EuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(TrimmedEuclideanDistancePointSetToPointSetMetricv4);

  /** Standard class type aliases. */
  using Self = TrimmedEuclideanDistancePointSetToPointSetMetricv4;
  using Superclass = EuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>;
  using SuperclassPointer = typename Superclass::Pointer;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TrimmedEuclideanDistancePointSetToPointSetMetricv4, EuclideanDistancePointSetToPointSetMetricv4 );

  /**  Type of the measure. */
  using MeasureType = typename Superclass::MeasureType;

  /**  Type of the parameters. */
  using ParametersType = typename Superclass::ParametersType;
  using ParametersValueType = typename Superclass::ParametersValueType;
  using NumberOfParametersType = typename Superclass::NumberOfParametersType;

  /**  Type of the derivative. */
  using DerivativeType = typename Superclass::DerivativeType;

  /** Transform types from Superclass*/
  using FixedTransformType = typename Superclass::FixedTransformType;
  using FixedTransformPointer = typename FixedTransformType::Pointer;
  using FixedTransformConstPointer = typename FixedTransformType::ConstPointer;
  using FixedInputPointType = typename Superclass::FixedInputPointType;
  using FixedOutputPointType = typename Superclass::FixedOutputPointType;
  using FixedTransformParametersType = typename Superclass::FixedTransformParametersType;

  using MovingTransformType = typename Superclass::MovingTransformType;
  using MovingTransformPointer = typename MovingTransformType::Pointer;
  using MovingTransformConstPointer = typename MovingTransformType::ConstPointer;
  using MovingInputPointType = typename Superclass::MovingInputPointType;
  using MovingOutputPointType = typename Superclass::MovingOutputPointType;
  using MovingTransformParametersType = typename Superclass::MovingTransformParametersType;

  using JacobianType = typename Superclass::JacobianType;
  using FixedTransformJacobianType = typename Superclass::FixedTransformJacobianType;
  using MovingTransformJacobianType = typename Superclass::MovingTransformJacobianType;

  using DisplacementFieldTransformType = typename Superclass::MovingDisplacementFieldTransformType;

  using ObjectType = typename Superclass::ObjectType;

  /** Dimension type */
  using DimensionType = typename Superclass::DimensionType;

  /**  Type of the fixed point set. */
  using FixedPointSetType = TFixedPointSet;
  using FixedPointSetPointer = typename TFixedPointSet::Pointer;
  using FixedPointSetConstPointer = typename TFixedPointSet::ConstPointer;
  using FixedPointType = typename TFixedPointSet::PointType;
  using FixedPointIdentifier = typename TFixedPointSet::PointIdentifier;
  using FixedPixelType = typename TFixedPointSet::PixelType;
  using FixedPointsContainer = typename TFixedPointSet::PointsContainer;

  static constexpr DimensionType FixedPointDimension = Superclass::FixedDimension;

  /**  Type of the moving point set. */
  using MovingPointSetType = TMovingPointSet;
  using MovingPointType = typename TMovingPointSet::PointType;
  using MovingPixelType = typename TMovingPointSet::PixelType;
  using MovingPointsContainer = typename TMovingPointSet::PointsContainer;

  static constexpr DimensionType MovingPointDimension = Superclass::MovingDimension;

  /**
   * typedefs for the data types used in the point set metric calculations.
   * It is assumed that the constants of the fixed point set, such as the
   * point dimension, are the same for the "common space" in which the metric
   * calculation occurs.
   */
  static constexpr DimensionType PointDimension = Superclass::FixedDimension;

  using PointType = FixedPointType;
  using PixelType = FixedPixelType;
  using CoordRepType = typename PointType::CoordRepType;
  using PointsContainer = FixedPointsContainer;
  using PointsConstIterator = typename PointsContainer::ConstIterator;
  using PointIdentifier = typename PointsContainer::ElementIdentifier;

  /** Typedef for points locator class to speed up finding neighboring points */
  using PointsLocatorType = PointsLocator< PointsContainer>;
  using NeighborsIdentifierType = typename PointsLocatorType::NeighborsIdentifierType;

  using FixedTransformedPointSetType = PointSet<FixedPixelType, Self::PointDimension >;
  using FixedTransformedPointsContainer = typename FixedTransformedPointSetType::PointsContainer;
  using FixedTransformedPointSetPointer = typename FixedTransformedPointSetType::Pointer;
  using FixedTransformedPointSetConstPointer = typename FixedTransformedPointSetType::ConstPointer;
  using MovingTransformedPointSetType = PointSet<MovingPixelType, Self::PointDimension >;

  using DerivativeValueType = typename DerivativeType::ValueType;
  using LocalDerivativeType = FixedArray<DerivativeValueType, Self::PointDimension >;

  /** Types for the virtual domain */
  using VirtualImageType = typename Superclass::VirtualImageType;
  using VirtualImagePointer = typename Superclass::VirtualImagePointer;
  using VirtualPixelType = typename Superclass::VirtualPixelType;
  using VirtualRegionType = typename Superclass::VirtualRegionType;
  using VirtualSizeType = typename Superclass::VirtualSizeType;
  using VirtualSpacingType = typename Superclass::VirtualSpacingType;
  using VirtualOriginType = typename Superclass::VirtualPointType;
  using VirtualPointType = typename Superclass::VirtualPointType;
  using VirtualDirectionType = typename Superclass::VirtualDirectionType;
  using VirtualRadiusType = typename Superclass::VirtualSizeType;
  using VirtualIndexType = typename Superclass::VirtualIndexType;
  using VirtualPointSetType = typename Superclass::VirtualPointSetType;
  using VirtualPointSetPointer = typename VirtualPointSetType::Pointer;
  using VirtualPointSetConstPointer = typename VirtualPointSetType::ConstPointer;
  using VirtualPointsContainer = typename VirtualPointSetType::PointsContainer;

  /**
   * This method returns the value of the metric based on the current
   * transformation(s).  This function can be redefined in derived classes
   * but many point set metrics follow the same structure---one iterates
   * through the points and, for each point a metric value is calculated.
   * The summation of these individual point metric values gives the total
   * value of the metric.  Note that this might not be applicable to all
   * point set metrics.  For those cases, the developer will have to redefine
   * the GetValue() function.
   */
  MeasureType GetValue() const override;

  /**
   * This method returns the derivative based on the current
   * transformation(s).  This function can be redefined in derived classes
   * but many point set metrics follow the same structure---one iterates
   * through the points and, for each point a derivative is calculated.
   * The set of all these local derivatives constitutes the total derivative.
   * Note that this might not be applicable to all point set metrics.  For
   * those cases, the developer will have to redefine the GetDerivative()
   * function.
   */
  void GetDerivative( DerivativeType & ) const override;

  /**
   * This method returns the derivative and value based on the current
   * transformation(s).  This function can be redefined in derived classes
   * but many point set metrics follow the same structure---one iterates
   * through the points and, for each point a derivative and value is calculated.
   * The set of all these local derivatives/values constitutes the total
   * derivative and value.  Note that this might not be applicable to all
   * point set metrics.  For those cases, the developer will have to redefine
   * the GetValue() and GetDerivative() functions.
   */
  void GetValueAndDerivative( MeasureType &, DerivativeType & ) const override;


  /**
   * Get/Set distance cutoff
   */
  itkGetMacro( DistanceCutoff, TInternalComputationValueType );
  itkSetMacro( DistanceCutoff, TInternalComputationValueType );

  /**
  * Get/Set the cut off percentile cut off value.
  */
  void SetPercentile( unsigned int percentile )
    {
    if( percentile > 0 && percentile <= 100 )
      {
      m_Percentile = percentile;
      }
    else
      {
      itkExceptionMacro( "Percentile value must belong to (0;100]." )
      }
    }
  itkGetMacro( Percentile, unsigned int );

  /**
  * Get/Set the cut off percentile cut off value.
  */
  void SetSamplingRate( double rate )
    {
    if( rate > 0 && rate <= 1 )
      {
      m_SamplingRate = rate;
      }
    else
      {
      itkExceptionMacro( "Sampling percentage value must belong to (0;1]." )
      }
    }
  itkGetMacro( SamplingRate, unsigned int );


protected:
  TrimmedEuclideanDistancePointSetToPointSetMetricv4();
  ~TrimmedEuclideanDistancePointSetToPointSetMetricv4() override = default;
  void PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Helper method allows for code reuse while skipping the metric value
   * calculation when appropriate */
  void MyCalculateValueAndDerivative( MeasureType & value, DerivativeType & derivative, bool calculateValue ) const;


private:
  /**
  * Cut off value to filter the number of points used to drive the registration
  * at each iteration. Points from the fixed point set whose error measurement
  * are within the given percentile are kept (i.e. if percentile is set to 80,
  * all the points below the 80th percentile are kept) to compute the error
  * and the derivative. Other points are ignored.
  */
  unsigned int m_Percentile;

  /**
   * Dsitance based cut off value, remove points larger than distance cutoff.
   * Can be used in conjunction with percentile filtering
   * */
  TInternalComputationValueType m_DistanceCutoff;

  /**
   * Use subsampling to compute gradient
   */
  double m_SamplingRate;

  //Create ranges over the point set for multithreaded computation of value and derivatives
  using PointIdentifierPair = std::pair<PointIdentifier, PointIdentifier>;
  using PointIdentifierRanges = std::vector<PointIdentifierPair>;
  const PointIdentifierRanges CreateRanges() const;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrimmedEuclideanDistancePointSetToPointSetMetricv4.hxx"
#endif

#endif
