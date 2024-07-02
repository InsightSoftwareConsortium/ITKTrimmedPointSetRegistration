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
#ifndef itkWeightedEuclideanDistancePointSetToPointSetMetricv4_h
#define itkWeightedEuclideanDistancePointSetToPointSetMetricv4_h

#include "itkPointSetToPointSetMetricv4.h"

namespace itk
{
/** \class WeightedEuclideanDistancePointSetToPointSetMetricv4
 * \brief Computes the Euclidan distance metric between two point sets.
 *
 * \ingroup TrimmedPointSetRegistration
 */
template<typename TFixedPointSet, typename TMovingPointSet = TFixedPointSet,
  class TInternalComputationValueType = double>
class ITK_TEMPLATE_EXPORT WeightedEuclideanDistancePointSetToPointSetMetricv4:
  public PointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(WeightedEuclideanDistancePointSetToPointSetMetricv4);

  /** Standard class type aliases. */
  using Self = WeightedEuclideanDistancePointSetToPointSetMetricv4;
  using Superclass = PointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet,
    TInternalComputationValueType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( WeightedEuclideanDistancePointSetToPointSetMetricv4, PointSetToPointSetMetricv4 );

  /** Types transferred from the base class */
  using MeasureType = typename Superclass::MeasureType;
  using DerivativeType = typename Superclass::DerivativeType;
  using LocalDerivativeType = typename Superclass::LocalDerivativeType;
  using PointType = typename Superclass::PointType;
  using PixelType = typename Superclass::PixelType;
  using PointIdentifier = typename Superclass::PointIdentifier;

  /**
   * Calculates the local metric value for a single point.
   */
  MeasureType GetLocalNeighborhoodValue( const PointType &, const PixelType & pixel = 0 ) const override;

  /**
   * Calculates the local value and derivative for a single point.
   */
  void GetLocalNeighborhoodValueAndDerivative( const PointType &,
    MeasureType &, LocalDerivativeType &, const PixelType & pixel = 0 ) const override;

  itkSetMacro(Weight, double);
  itkGetMacro(Weight, double);
protected:
  WeightedEuclideanDistancePointSetToPointSetMetricv4() = default;
  ~WeightedEuclideanDistancePointSetToPointSetMetricv4() override = default;

  /** PrintSelf function */
  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  double ComputeWeight(const double distance) const
    {
    return exp(-(distance*distance) / (m_Weight*m_Weight) );
    }
  double m_Weight;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWeightedEuclideanDistancePointSetToPointSetMetricv4.hxx"
#endif

#endif
