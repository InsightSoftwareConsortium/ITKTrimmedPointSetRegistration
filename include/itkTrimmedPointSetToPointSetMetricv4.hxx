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
#ifndef itkTrimmedPointSetToPointSetMetric_hxx
#define itkTrimmedPointSetToPointSetMetric_hxx

#include "itkTrimmedPointSetToPointSetMetricv4.h"
#include "itkIdentityTransform.h"
#include "itkCompensatedSummation.h"

namespace itk
{

/** Constructor */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::TrimmedPointSetToPointSetMetricv4()
{
  this->m_Metric = nullptr;    // has to be provided by the user.
  this->m_Percentile = 100;
  this->m_DistanceCutoff = NumericTraits<TInternalComputationValueType>::max();
}

/** Initialize the metric */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::Initialize()
{
  if ( !this->m_Metric )
    {
    itkExceptionMacro( "Decorated Metric is not present" );
    }
  m_Metric->Initialize();
}


template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::MeasureType
TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValue() const
{
  m_Metric->InitializeForIteration();

  FixedPointSetPointer fixedPointSet = m_Metric->GetFixedPointSet();
  FixedTransformedPointSetPointer fixedTransformedPointSet = m_Metric->GetFixedTransformedPointSet();
  VirtualPointSetPointer virtualTransformedPointSet = m_Metric->GetVirtualTransformedPointSet();

  PointsConstIterator It = fixedTransformedPointSet->GetPoints()->Begin();

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  if( virtualTransformedPointSet->GetNumberOfPoints() != fixedTransformedPointSet->GetNumberOfPoints() )
    {
    itkExceptionMacro("Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet.");
    }
  PointsConstIterator virtualIt = virtualTransformedPointSet->GetPoints()->Begin();
  PointIdentifier numberOfFixedPoints = fixedTransformedPointSet->GetNumberOfPoints();
  std::vector<MeasureType> distances(numberOfFixedPoints, NumericTraits<MeasureType>::max());
  size_t index = 0;
  while( It != fixedTransformedPointSet->GetPoints()->End() )
    {
    /* Verify the virtual point is in the virtual domain.
     * If user hasn't defined a virtual space, and the active transform is not
     * a displacement field transform type, then this will always return true. */
    if( ! m_Metric->IsInsideVirtualDomain( virtualIt.Value() ) )
      {
      ++It;
      ++virtualIt;
      continue;
      }

    PixelType pixel;
    NumericTraits<PixelType>::SetLength( pixel, 1 );
    if( m_Metric->m_UsePointSetData )
      {
      bool doesPointDataExist = fixedPointSet->GetPointData( It.Index(), &pixel );
      if( ! doesPointDataExist )
        {
        itkExceptionMacro( "The corresponding data for point " << It.Value() << " (pointId = " << It.Index() << ") does not exist." );
        }
      }
    distances[index] = m_Metric->GetLocalNeighborhoodValue( It.Value(), pixel );
    ++virtualIt;
    ++It;
    ++index;
    }
  DerivativeType derivative;
  // `valueSum` default value is set to max in `VerifyNumberOfValidPoints`
  // if there is no valid point.
  MeasureType valueSum = 0.0;
  if( this->VerifyNumberOfValidPoints( valueSum, derivative ) )
    {
    size_t last_index = m_Metric->GetNumberOfValidPoints();
    if( m_Percentile < 100)
      {
      std::sort(distances.begin(), distances.end());
      last_index = this->m_NumberOfValidPoints * m_Percentile / 100;
      }
    if( last_index < 1 )
      {
      itkExceptionMacro( "Percentile too small. No valid point in selected percentile." );
      }

    unsigned int nDistances = 0;
    for( typename std::vector<MeasureType>::iterator it = distances.begin();
         it <= distances.begin() + last_index; it++ )
      {
      if( (*it) < m_DistanceCutoff )
        {
        valueSum += (*it);
        ++nDistances;
        }
      }
    if( nDistances < 1 )
      {
      itkExceptionMacro( "DistanceCutoff too small. No valid point in selected percentile." );
      }
    valueSum /= nDistances;
    }
  this->m_Value = valueSum;

  return valueSum;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetDerivative( DerivativeType & derivative ) const
{
  MeasureType value = NumericTraits<MeasureType>::ZeroValue();
  this->CalculateValueAndDerivative( value, derivative, false );
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValueAndDerivative( MeasureType & value, DerivativeType & derivative ) const
{
  this->CalculateValueAndDerivative( value, derivative, true );
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::CalculateValueAndDerivative( MeasureType & calculatedValue, DerivativeType & derivative, bool calculateValue ) const
{
  m_Metric->InitializeForIteration();

  FixedPointSetPointer fixedPointSet = m_Metric->GetFixedPointSet();
  FixedTransformedPointSetPointer fixedTransformedPointSet = m_Metric->GetFixedTransformedPointSet();
  VirtualPointSetPointer virtualTransformedPointSet = m_Metric->GetVirtualTransformedPointSet();


  derivative.SetSize( m_Metric->GetNumberOfParameters() );
  if( ! this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
    {
    derivative.SetSize( PointDimension * fixedTransformedPointSet->GetNumberOfPoints() );
    }
  derivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

  CompensatedSummation<MeasureType> value;
  MovingTransformJacobianType  jacobian( MovingPointDimension, m_Metric->GetNumberOfLocalParameters() );
  MovingTransformJacobianType  jacobianCache;

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  PointIdentifier numberOfFixedPoints = fixedTransformedPointSet->GetNumberOfPoints();
  if( virtualTransformedPointSet->GetNumberOfPoints() != numberOfFixedPoints )
    {
    itkExceptionMacro( "Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet." );
    }
  PointsConstIterator virtualIt = virtualTransformedPointSet->GetPoints()->Begin();
  PointsConstIterator It = fixedTransformedPointSet->GetPoints()->Begin();
  PointsConstIterator end = fixedTransformedPointSet->GetPoints()->End();
  using ValueType = std::pair<MeasureType, DerivativeType>;
  std::vector<ValueType> values(numberOfFixedPoints);
  size_t index = 0;
  while( It != end )
    {
    DerivativeType localTransformDerivative( m_Metric->GetNumberOfLocalParameters() );
    localTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );
    MeasureType pointValue = NumericTraits<MeasureType>::ZeroValue();
    LocalDerivativeType pointDerivative;

    /* Verify the virtual point is in the virtual domain.
     * If user hasn't defined a virtual space, and the active transform is not
     * a displacement field transform type, then this will always return true. */
    if( ! m_Metric->IsInsideVirtualDomain( virtualIt.Value() ) )
      {
      ++It;
      ++virtualIt;
      values[index].first = NumericTraits<MeasureType>::max();
      continue;
      }

    PixelType pixel;
    NumericTraits<PixelType>::SetLength( pixel, 1 );
    if( m_Metric->m_UsePointSetData )
      {
      bool doesPointDataExist = fixedPointSet->GetPointData( It.Index(), &pixel );
      if( ! doesPointDataExist )
        {
        itkExceptionMacro( "The corresponding data for point " << It.Value() << " (pointId = " << It.Index() << ") does not exist." );
        }
      }

    if( calculateValue )
      {
      m_Metric->GetLocalNeighborhoodValueAndDerivative( It.Value(), pointValue, pointDerivative, pixel );
      values[index].first = pointValue;
      }
    else
      {
      pointDerivative = m_Metric->GetLocalNeighborhoodDerivative( It.Value(), pixel );
      }

    if( m_Metric->GetCalculateValueAndDerivativeInTangentSpace() )
      {
      jacobian.Fill( 0.0 );
      for( DimensionType d = 0; d < MovingPointDimension; d++ )
        {
        jacobian(d, d) = 1.0;
        }
      }
    else
      {
      m_Metric->GetMovingTransform()->
        ComputeJacobianWithRespectToParametersCachedTemporaries( virtualIt.Value(),
                                                                 jacobian,
                                                                 jacobianCache );
      }

    for( NumberOfParametersType par = 0; par < m_Metric->GetNumberOfLocalParameters(); par++ )
      {
      for( DimensionType d = 0; d < PointDimension; ++d )
        {
        localTransformDerivative[par] += jacobian(d, par) * pointDerivative[d];
        }
      }

    // For local-support transforms, store the per-point result
    if( m_Metric->HasLocalSupport() || m_Metric->m_CalculateValueAndDerivativeInTangentSpace )
      {
      if( m_Metric->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
        {
        m_Metric->StorePointDerivative( virtualIt.Value(), localTransformDerivative, derivative );
        }
      else
        {
        for( NumberOfParametersType par = 0; par < m_Metric->GetNumberOfLocalParameters(); par++ )
          {
          derivative[m_Metric->GetNumberOfLocalParameters() * It.Index() + par] = localTransformDerivative[par];
          }
        }
      }
    else
      {
      values[index].second = localTransformDerivative;
      }

    ++It;
    ++virtualIt;
    ++index;
    }

  // `valueSum` default value is set to max in `VerifyNumberOfValidPoints`
  // if there is no valid point.
  MeasureType valueSum = 0.0;
  if( this->VerifyNumberOfValidPoints( valueSum, derivative ) )
    {
    size_t last_index = m_Metric->GetNumberOfValidPoints();
    if( m_Percentile < 100)
      {
      std::sort( values.begin(), values.end(), [](ValueType a, ValueType b)
        { return a.first < b.first ? true : false; });

      last_index = m_Metric->GetNumberOfValidPoints() * m_Percentile / 100;
      }
    if( last_index < 1 )
      {
      itkExceptionMacro( "Percentile too small. No valid point in selected percentile." );
      }

    unsigned int nDistances = 0;
    for( typename std::vector<ValueType>::iterator it = values.begin();
         it <= values.begin() + last_index; it++ )
      {
      ValueType &el = *it;
      if( el.first < m_DistanceCutoff )
        {
        valueSum += el.first;
        ++nDistances;
        if( ! m_Metric->HasLocalSupport() && ! m_Metric->GetCalculateValueAndDerivativeInTangentSpace() )
          {
          for( NumberOfParametersType par = 0; par < m_Metric->GetNumberOfLocalParameters(); par++ )
            {
            derivative[par] += el.second[par];
            }
          }
        }
      }
    if( nDistances < 1 )
      {
      itkExceptionMacro( "DistanceCutoff too small. No valid point in selected percentile." );
      }

    valueSum /= nDistances;
    derivative /= nDistances;
    }

  calculatedValue = valueSum;
  m_Metric->m_Value = valueSum;
  std::cout<<"Get value and derivative"<<std::endl;
}



/** PrintSelf */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedPointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Trimmed Metric: " << this->m_Metric.GetPointer() << std::endl;
  m_Metric->PrintSelf();
}
} // end namespace itk

#endif
