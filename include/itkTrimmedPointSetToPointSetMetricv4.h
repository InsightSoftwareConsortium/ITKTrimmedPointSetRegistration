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
#ifndef itkTrimmedPointSetToPointSetMetricv4_h
#define itkTrimmedPointSetToPointSetMetricv4_h

#include "itkPointSetToPointSetMetricv4.h"

namespace itk
{
/** \class TrimmedPointSetToPointSetMetricv4
 * \brief Computes similarity between two point sets.
 *
 * \ingroup ITKMetricsv4
 */

template<typename TFixedPointSet,  typename TMovingPointSet,
  class TInternalComputationValueType = double>
class ITK_TEMPLATE_EXPORT TrimmedPointSetToPointSetMetricv4
: public PointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(TrimmedPointSetToPointSetMetricv4);

  /** Standard class type aliases. */
  using Self = TrimmedPointSetToPointSetMetricv4;
  using Superclass = PointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

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
  using FixedTransformPointer = typename Superclass::FixedTransformPointer;
  using FixedInputPointType = typename Superclass::FixedInputPointType;
  using FixedOutputPointType = typename Superclass::FixedOutputPointType;
  using FixedTransformParametersType = typename Superclass::FixedTransformParametersType;

  using MovingTransformType = typename Superclass::MovingTransformType;
  using MovingTransformPointer = typename Superclass::MovingTransformPointer;
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
  using VirtualPointSetPointer = typename Superclass::VirtualPointSetPointer;

  /** Set fixed point set*/
  void SetFixedObject( const ObjectType *object ) override
    {
    m_Metric.SetFixedObject();
    }

  /** Set moving point set*/
  void SetMovingObject( const ObjectType *object ) override
    {
    m_Metric.SetMovingObject();
    }

  /** Get/Set the fixed pointset.  */
  itkSetConstObjectMacro( FixedPointSet, FixedPointSetType );
  itkGetConstObjectMacro( FixedPointSet, FixedPointSetType );

  /** Get the fixed transformed point set.  */
  itkGetModifiableObjectMacro( FixedTransformedPointSet, FixedTransformedPointSetType );

  /** Get/Set the moving point set.  */
  itkSetConstObjectMacro( MovingPointSet, MovingPointSetType );
  itkGetConstObjectMacro( MovingPointSet, MovingPointSetType );

  /** Get the moving transformed point set.  */
  itkGetModifiableObjectMacro( MovingTransformedPointSet, MovingTransformedPointSetType );

  /**
   * For now return the number of points used in the value/derivative calculations.
   */
  SizeValueType GetNumberOfComponents() const;

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
  virtual MeasureType GetLocalNeighborhoodValue( const PointType & p, const PixelType & pixel ) const
    {
    return m_Metric.GetLocalNeighborhoodValue(p, pixel);
    }

  /**
   * Calculates the local derivative for a single point. The \c PixelType may or
   * may not be used.  See class description for further explanation.
   */
  virtual LocalDerivativeType GetLocalNeighborhoodDerivative( const PointType &p, const PixelType & pixel ) const
    {
    return m_Metric.GetLocalNeighborhoodDerivative(p, pixel);
    }

  /**
   * Calculates the local value/derivative for a single point.  The \c PixelType may or
   * may not be used.  See class description for further explanation.
   */
  virtual void GetLocalNeighborhoodValueAndDerivative( const PointType &p,
    MeasureType &m, LocalDerivativeType &l, const PixelType & pixel ) const
    {
    m_Metric.GetLocalNeighborhoodDerivative(p, m, l, pixel);
    }

  /**
   * Get the virtual point set, derived from the fixed point set.
   * If the virtual point set has not yet been derived, it will be
   * in this call. */
  const VirtualPointSetType * GetVirtualTransformedPointSet() const
    {
    return m_Metric.GetVirtualTransformedPointSet();
    }

  /**
   * Initialize the metric by making sure that all the components
   *  are present and plugged together correctly.
   */
  void Initialize() override
    {
    m_Metric.Initialize();
    }

  bool SupportsArbitraryVirtualDomainSamples() const override
    {
    return m_Metric.SupportsArbitraryVirtualDomainSamples();
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
  itkSetMacro( StoreDerivativeAsSparseFieldForLocalSupportTransforms, bool );
  itkGetConstMacro( StoreDerivativeAsSparseFieldForLocalSupportTransforms, bool );
  itkBooleanMacro( StoreDerivativeAsSparseFieldForLocalSupportTransforms );

  /**
   *
   */
  itkSetMacro( CalculateValueAndDerivativeInTangentSpace, bool );
  itkGetConstMacro( CalculateValueAndDerivativeInTangentSpace, bool );
  itkBooleanMacro( CalculateValueAndDerivativeInTangentSpace );

protected:
  TrimmedPointSetToPointSetMetricv4();
  ~TrimmedPointSetToPointSetMetricv4() override = default;
  void PrintSelf( std::ostream & os, Indent indent ) const override;

  typename FixedPointSetType::ConstPointer                m_FixedPointSet;
  mutable typename FixedTransformedPointSetType::Pointer  m_FixedTransformedPointSet;

  mutable typename PointsLocatorType::Pointer             m_FixedTransformedPointsLocator;

  typename MovingPointSetType::ConstPointer               m_MovingPointSet;
  mutable typename MovingTransformedPointSetType::Pointer m_MovingTransformedPointSet;

  mutable typename PointsLocatorType::Pointer             m_MovingTransformedPointsLocator;

  /** Holds the fixed points after transformation into virtual domain. */
  mutable VirtualPointSetPointer                          m_VirtualTransformedPointSet;


  /**
   * Prepare point sets for use. */
  virtual void InitializePointSets() const
    {
    m_Metric.InitializePointSets();
    }

  /**
   * Initialize to prepare for a particular iteration, generally
   * an iteration of optimization. Distinct from Initialize()
   * which is a one-time initialization. */
  virtual void InitializeForIteration() const
    {
    m_Metric.InitializeForIteration();
    }

  /**
   * Determine the number of valid fixed points. A fixed point
   * is valid if, when transformed into the virtual domain using
   * the inverse of the FixedTransform, it is within the defined
   * virtual domain bounds. */
  virtual SizeValueType CalculateNumberOfValidFixedPoints() const
    {
    m_Metric.CalculateNumberOfValidFixedPoints();
    }

  /** Helper method allows for code reuse while skipping the metric value
   * calculation when appropriate */
  void CalculateValueAndDerivative( MeasureType & value, DerivativeType & derivative, bool calculateValue ) const;

  /**
   * Warp the fixed point set into the moving domain based on the fixed transform,
   * passing through the virtual domain and storing a virtual domain set.
   * Note that the warped moving point set is of type FixedPointSetType since the transform
   * takes the points from the fixed to the moving domain.
   */
  void TransformFixedAndCreateVirtualPointSet() const
    {
    m_Metric.TransformFixedAndCreateVirtualPointSet();
    }

  /**
   * Warp the moving point set based on the moving transform.  Note that the
   * warped moving point set is of type FixedPointSetType since the transform
   * takes the points from the moving to the fixed domain.
   * FIXME: needs update.
   */
  void TransformMovingPointSet() const
    {
    m_Metric.TransformMovingPointSet();
    }

  /**
   * Build point locators for the fixed and moving point sets to speed up
   * derivative and value calculations.
   */
  void InitializePointsLocators() const
    {
    m_Metric.InitializePointsLocators();
    }

  /**
   * Store a derivative from a single point in a field.
   * Only relevant when active transform has local support.
   */
  void StorePointDerivative( const VirtualPointType &, const DerivativeType &, DerivativeType & ) const;

  using MetricCategoryType = typename Superclass::MetricCategoryType;

  /** Get metric category */
  MetricCategoryType GetMetricCategory() const override
    {
    return m_Metric.GetMetricCategory();
    }


private:
  PointSetToPointSetMetricv4::Pointer m_Metric;
  mutable bool m_MovingTransformPointLocatorsNeedInitialization;
  mutable bool m_FixedTransformPointLocatorsNeedInitialization;

  // Flag to keep track of whether a warning has already been issued
  // regarding the number of valid points.
  mutable bool m_HaveWarnedAboutNumberOfValidPoints;

  // Flag to store derivatives at fixed point locations with the rest being zero gradient
  // (default = true).
  bool m_StoreDerivativeAsSparseFieldForLocalSupportTransforms;

  mutable ModifiedTimeType m_MovingTransformedPointSetTime;
  mutable ModifiedTimeType m_FixedTransformedPointSetTime;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTrimmedPointSetToPointSetMetricv4.hxx"
#endif

#endif
