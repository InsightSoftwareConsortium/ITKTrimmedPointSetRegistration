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
#ifndef itkTrimmedPointSetToPointSetMetricv4_h
#define itkTrimmedPointSetToPointSetMetricv4_h

#include "itkEuclideanDistancePointSetToPointSetMetricv4.h"

namespace itk
{
/** \class TrimmedPointSetToPointSetMetricv4
 * \brief Computes similarity between two point sets.
 *
 * \ingroup TrimmedPointSetRegistration
 */

template<typename TFixedPointSet,  typename TMovingPointSet,
  class TInternalComputationValueType = double>
class ITK_TEMPLATE_EXPORT TrimmedPointSetToPointSetMetricv4
: public EuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(TrimmedPointSetToPointSetMetricv4);

  /** Standard class type aliases. */
  using Self = TrimmedPointSetToPointSetMetricv4;
  using Superclass = EuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>;
  using SuperclassPointer = typename Superclass::Pointer;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TrimmedPointSetToPointSetMetricv4, PointSetToPointSetMetricv4 );

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

  /** Set fixed point set*/
  void SetFixedObject( const ObjectType *object ) override
    {
    Superclass::SetFixedObject(object);
    m_Metric->SetFixedObject(object);
    }

  /** Set moving point set*/
  void SetMovingObject( const ObjectType *object ) override
    {
    Superclass::SetMovingObject(object);
    m_Metric->SetMovingObject(object);
    }

  /** Get/Set the fixed pointset.  */
  void SetFixedPointSet(const FixedPointSetType *set) override
    {
    Superclass::SetFixedPointSet(set);
    m_Metric->SetFixedPointSet(set);
    }

  const FixedPointSetType * GetFixedPointSet() const override
    {
    return m_Metric->GetFixedPointSet();
    }


  /** Get the fixed transformed point set.  */
  FixedTransformedPointSetType * GetFixedTransformedPointSet() override
    {
    return m_Metric->GetFixedTransformedPointSet();
    }

  /** Get/Set the moving point set.  */
  void SetMovingPointSet(const MovingPointSetType *set) override
    {
    Superclass::SetMovingPointSet(set);
    m_Metric->SetMovingPointSet(set);
    }

  const MovingPointSetType * GetMovingPointSet() const override
    {
    return m_Metric->GetMovingPointSet();
    }


  /** Get the moving transformed point set.  */
  MovingTransformedPointSetType * GetMovingTransformedPointSet() override
    {
    return m_Metric->GetFixedTransformedPointSet();
    }

  /**
   * For now return the number of points used in the value/derivative calculations.
   */
  SizeValueType GetNumberOfComponents() const
    {
    return m_Metric->GetNumberOfComponents();
    }


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
   * Function to be defined in the appropriate derived classes.  Calculates
   * the local metric value for a single point.  The \c PixelType may or
   * may not be used.  See class description for further explanation.
   */
  MeasureType GetLocalNeighborhoodValue( const PointType & p, const PixelType & pixel ) const override
    {
    return m_Metric->GetLocalNeighborhoodValue(p, pixel);
    }

  /**
   * Calculates the local derivative for a single point. The \c PixelType may or
   * may not be used.  See class description for further explanation.
   */
  LocalDerivativeType GetLocalNeighborhoodDerivative( const PointType &p, const PixelType & pixel ) const override
    {
    return m_Metric->GetLocalNeighborhoodDerivative(p, pixel);
    }

  /**
   * Calculates the local value/derivative for a single point.  The \c PixelType may or
   * may not be used.  See class description for further explanation.
   */
  void GetLocalNeighborhoodValueAndDerivative( const PointType &p,
    MeasureType &m, LocalDerivativeType &l, const PixelType & pixel ) const override
    {
    m_Metric->GetLocalNeighborhoodValueAndDerivative(p, m, l, pixel);
    }

  /**
   * Get the virtual point set, derived from the fixed point set.
   * If the virtual point set has not yet been derived, it will be
   * in this call. */
  const VirtualPointSetType * GetVirtualTransformedPointSet() const
    {
    return m_Metric->GetVirtualTransformedPointSet();
    }

  /**
   * Initialize the metric by making sure that all the components
   *  are present and plugged together correctly.
   */
  void Initialize() override;

  bool SupportsArbitraryVirtualDomainSamples() const override
    {
    return m_Metric->SupportsArbitraryVirtualDomainSamples();
    }

  /**
   * By default, the point set metric derivative for a displacement field transform
   * is stored by saving the gradient for every voxel in the displacement field (see
   * the function StorePointDerivative()).  Since the "fixed points" will typically
   * constitute a sparse set, this means that the field will have zero gradient values
   * at every voxel that doesn't have a corresponding point.  This might cause additional
   * computation time for certain transforms (e.g. B-spline SyN). To avoid this, this
   * option permits storing the point derivative only at the fixed point locations.
   * If this variable is set to false, then the derivative array will be of length
   * = PointDimension * m_FixedPointSet->GetNumberOfPoints().
   */
  void SetStoreDerivativeAsSparseFieldForLocalSupportTransforms(bool sparse) override
    {
    Superclass::SetStoreDerivativeAsSparseFieldForLocalSupportTransforms(sparse);
    m_Metric->SetStoreDerivativeAsSparseFieldForLocalSupportTransforms(sparse);
    }
  bool GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() const override
    {
    return m_Metric->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms();
    }
  void StoreDerivativeAsSparseFieldForLocalSupportTransformsOn() override
    {
    Superclass::StoreDerivativeAsSparseFieldForLocalSupportTransformsOn();
    m_Metric->StoreDerivativeAsSparseFieldForLocalSupportTransformsOn();
    }
  void StoreDerivativeAsSparseFieldForLocalSupportTransformsOff() override
    {
    Superclass::StoreDerivativeAsSparseFieldForLocalSupportTransformsOff();
    m_Metric->StoreDerivativeAsSparseFieldForLocalSupportTransformsOff();
    }

  /**
   *
   */
  void SetCalculateValueAndDerivativeInTangentSpace(bool sparse) override
    {
    Superclass::SetCalculateValueAndDerivativeInTangentSpace(sparse);
    m_Metric->SetCalculateValueAndDerivativeInTangentSpace(sparse);
    }
  bool GetCalculateValueAndDerivativeInTangentSpace() const override
    {
    return m_Metric->GetCalculateValueAndDerivativeInTangentSpace();
    }
  void CalculateValueAndDerivativeInTangentSpaceOn() override
    {
    Superclass::CalculateValueAndDerivativeInTangentSpaceOn();
    m_Metric->CalculateValueAndDerivativeInTangentSpaceOn();
    }
  void CalculateValueAndDerivativeInTangentSpaceOff() override
    {
    Superclass::CalculateValueAndDerivativeInTangentSpaceOff();
    m_Metric->CalculateValueAndDerivativeInTangentSpaceOff();
    }

  /*
  void SetMetric(Superclass::Pointer metric)
    {
      m_Metric = metric;
      this->SetFixedPointSet( metric->GetFixedPointSet() );
      this->SetFixedTransform( metric->GetFixedTransform() );
      this->SetMovingPointSet( metric->GetMovingPointSet() );
      this->SetMovingTransform( metric->GetMovingTransform() );
      this->SetStoreDerivativeAsSparseFieldForLocalSupportTransforms(
                      metric->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() );
      this->SetCalculateValueAndDerivativeInTangentSpace(
                      metric->GetCalculateValueAndDerivativeInTangentSpace() );
    }
  */
  itkSetObjectMacro(Metric, Superclass);
  itkGetConstObjectMacro(Metric, Superclass)


  void SetFixedTransform( FixedTransformType *fixed) override
    {
    Superclass::SetFixedTransform(fixed);
    m_Metric->SetFixedTransform(fixed);
    }

  void SetMovingTransform(MovingTransformType *moving) override
    {
    Superclass::SetMovingTransform(moving);
    m_Metric->SetMovingTransform(moving);
    }

  /**
   * Get/Set distance cutoff
   */
  itkGetMacro( DistanceCutoff, TInternalComputationValueType );

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


protected:
  TrimmedPointSetToPointSetMetricv4();
  ~TrimmedPointSetToPointSetMetricv4() override = default;
  void PrintSelf( std::ostream & os, Indent indent ) const override;


  /** Helper method allows for code reuse while skipping the metric value
   * calculation when appropriate */
  void CalculateValueAndDerivative( MeasureType & value, DerivativeType & derivative, bool calculateValue ) const;

  /**
   * Prepare point sets for use. */
  void InitializePointSets() const override
    {
    Superclass::InitializePointSets();
    //m_Metric->InitializePointSets();
    }

  /**
   * Initialize to prepare for a particular iteration, generally
   * an iteration of optimization. Distinct from Initialize()
   * which is a one-time initialization. */
  void InitializeForIteration() const override
    {
    Superclass::InitializeForIteration();
    //m_Metric->InitializeForIteration();
    }


private:
  //Decorated metric to be used
  SuperclassPointer m_Metric;

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
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrimmedPointSetToPointSetMetricv4.hxx"
#endif

#endif
