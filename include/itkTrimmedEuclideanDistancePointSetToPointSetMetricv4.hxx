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
#ifndef itkTrimmedEuclideanDistancePointSetToPointSetMetricv4_hxx
#define itkTrimmedEuclideanDistancePointSetToPointSetMetricv4_hxx

#include "itkTrimmedEuclideanDistancePointSetToPointSetMetricv4.h"
#include "itkIdentityTransform.h"
#include "itkCompensatedSummation.h"

#include "itkTimeProbe.h"

namespace itk
{

/** Constructor */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::TrimmedEuclideanDistancePointSetToPointSetMetricv4()
{
  this->m_Percentile = 100;
  this->m_DistanceCutoff = NumericTraits<TInternalComputationValueType>::max();
}


template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::MeasureType
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValue() const
{
  using FixedTransformedVectorContainer = typename FixedPointsContainer::STLContainerType;
  using VirtualVectorContainer =  typename VirtualPointsContainer::STLContainerType;
  
  this->InitializeForIteration();

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  PointIdentifier numberOfFixedPoints = this->GetFixedTransformedPointSet()->GetNumberOfPoints();
  if( this->GetVirtualTransformedPointSet()->GetNumberOfPoints() != numberOfFixedPoints )
    {
    itkExceptionMacro("Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet.");
    }

  std::vector<MeasureType> distances(numberOfFixedPoints, NumericTraits<MeasureType>::max());

  {
  const VirtualVectorContainer &virtualTransformedPointSet = 
    this->GetVirtualTransformedPointSet()->GetPoints()->CastToSTLConstContainer();  
  const FixedTransformedVectorContainer &fixedTransformedPointSet = 
    this->GetFixedTransformedPointSet()->GetPoints()->CastToSTLConstContainer();
  std::function< void(FixedPointIdentifier) > collectNeighborhoodValues =
          [&distances, this, &virtualTransformedPointSet, &fixedTransformedPointSet](FixedPointIdentifier index)
    {
    FixedPointType virtualTransformedPoint = virtualTransformedPointSet[index];
    FixedPointType fixedTransformedPoint =  fixedTransformedPointSet[index];
    if( this->IsInsideVirtualDomain( virtualTransformedPoint ) )
      {
      PixelType pixel;
      NumericTraits<PixelType>::SetLength( pixel, 1 );
      if( this->m_UsePointSetData )
        {
        bool doesPointDataExist = this->GetFixedPointSet()->GetPointData( index, &pixel );
        if( ! doesPointDataExist )
          {
          itkExceptionMacro( "The corresponding data for point (pointId = " << index << ") does not exist." );
          }
      }

      distances[ index ] = this->GetLocalNeighborhoodValue( fixedTransformedPoint, pixel );

      }

    };
  MultiThreaderBase::New()->ParallelizeArray( (FixedPointIdentifier) 0, 
                          (FixedPointIdentifier) fixedTransformedPointSet.size(), 
                           collectNeighborhoodValues, nullptr );
  }

  DerivativeType derivative;
  // `valueSum` default value is set to max in `VerifyNumberOfValidPoints`
  // if there is no valid point.
  MeasureType valueSum = 0.0;
  if( this->VerifyNumberOfValidPoints( valueSum, derivative ) )
    {
    size_t last_index = this->GetNumberOfValidPoints();
    if( m_Percentile < 100)
      {
      std::sort(distances.begin(), distances.end());
      last_index = this->GetNumberOfValidPoints() * m_Percentile / 100;
      }
    if( last_index < 1 )
      {
      itkExceptionMacro( "Percentile too small. No valid point in selected percentile." );
      }

    unsigned int nDistances = 0;
    for( typename std::vector<MeasureType>::iterator it = distances.begin();
         it < distances.begin() + last_index; it++ )
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
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetDerivative( DerivativeType & derivative ) const
{
  MeasureType value = NumericTraits<MeasureType>::ZeroValue();
  this->MyCalculateValueAndDerivative( value, derivative, false );
}


template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValueAndDerivative( MeasureType & value, DerivativeType & derivative ) const
{
  this->MyCalculateValueAndDerivative( value, derivative, true );
}


template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::MyCalculateValueAndDerivative( MeasureType & calculatedValue, DerivativeType & derivative, bool calculateValue ) const
{
  /*
  itk::TimeProbe clock1;
  itk::TimeProbe clock2;
  */
  struct PointDerivativeStorage
    {
    FixedPointIdentifier index;
    DerivativeType derivative;
    };
  using ValueType = std::pair<MeasureType, PointDerivativeStorage>;
  
  using FixedTransformedVectorContainer = typename FixedPointsContainer::STLContainerType;
  using VirtualVectorContainer =  typename VirtualPointsContainer::STLContainerType;

  this->InitializeForIteration();
  
  derivative.SetSize( this->GetNumberOfParameters() );
  if( ! this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
    {
    derivative.SetSize( PointDimension * this->GetFixedTransformedPointSet()->GetNumberOfPoints() );
    }
  derivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

  CompensatedSummation<MeasureType> value;

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  if( this->GetVirtualTransformedPointSet()->GetNumberOfPoints() !=
      this->GetFixedTransformedPointSet()->GetNumberOfPoints() )
    {
    itkExceptionMacro( "Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet." );
    }
  std::vector<ValueType> values( this->GetFixedTransformedPointSet()->GetNumberOfPoints() );
  
  //clock1.Start();
  {
  const VirtualVectorContainer &virtualTransformedPointSet = 
    this->GetVirtualTransformedPointSet()->GetPoints()->CastToSTLConstContainer();  
  const FixedTransformedVectorContainer &fixedTransformedPointSet = 
    this->GetFixedTransformedPointSet()->GetPoints()->CastToSTLConstContainer();

  //Collect derviatives at each point
  std::function< void(FixedPointIdentifier) > collectNeighborhoodValues =
       [&values, this, &calculateValue, &derivative, &virtualTransformedPointSet, &fixedTransformedPointSet]
       (FixedPointIdentifier index)
    {
    /* Verify the virtual point is in the virtual domain.
     * If user hasn't defined a virtual space, and the active transform is not
     * a displacement field transform type, then this will always return true. */
    values[index].first = 0;
    values[index].second.index = index;
    values[index].second.derivative.SetSize( this->GetNumberOfLocalParameters() );
    values[index].second.derivative.Fill(0);
    if( !this->IsInsideVirtualDomain( virtualTransformedPointSet[index] ) )
      {
      return;
      }
    PixelType pixel;
    NumericTraits<PixelType>::SetLength( pixel, 1 );
    if( this->m_UsePointSetData )
      {
      bool doesPointDataExist = this->GetFixedTransformedPointSet()->GetPointData( index, &pixel );
      if( ! doesPointDataExist )
        {
        itkExceptionMacro( "The corresponding data for point " << index << ") does not exist." );
        }
      }
    LocalDerivativeType pointDerivative;
    if( calculateValue || m_Percentile < 100 || 
        m_DistanceCutoff < NumericTraits<TInternalComputationValueType>::max() )
      {
      MeasureType pointValue = NumericTraits<MeasureType>::ZeroValue();
      this->GetLocalNeighborhoodValueAndDerivative( 
          fixedTransformedPointSet[index], pointValue, pointDerivative, pixel );
      values[index].first = pointValue;
      }
    else
      {
      pointDerivative = this->GetLocalNeighborhoodDerivative( fixedTransformedPointSet[index], pixel );
      }
            
    if( !this->GetCalculateValueAndDerivativeInTangentSpace() )
      {
      thread_local MovingTransformJacobianType jacobian( MovingPointDimension, this->GetNumberOfLocalParameters() );
      thread_local MovingTransformJacobianType jacobianCache;
      this->GetMovingTransform()->ComputeJacobianWithRespectToParametersCachedTemporaries(  
          virtualTransformedPointSet[index], jacobian, jacobianCache );
      for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
        {
        for( DimensionType d = 0; d < PointDimension; ++d )
          {
          values[index].second.derivative[par] += jacobian(d, par) * pointDerivative[d];
          }
        }
      }
    else
      {
      for( DimensionType d = 0; d < PointDimension; ++d )
        {
        values[index].second.derivative[d] = pointDerivative[d];
        }
      }  
    };
    
  MultiThreaderBase::New()->ParallelizeArray( (FixedPointIdentifier) 0, 
                          (FixedPointIdentifier) fixedTransformedPointSet.size(), 
                           collectNeighborhoodValues, nullptr );
  }
  //clock1.Stop();
  
  //clock2.Start();
  // `valueSum` default value is set to max in `VerifyNumberOfValidPoints`
  // if there is no valid point.
  MeasureType initialValueSum = 0;
  calculatedValue = initialValueSum;
  this->m_Value = initialValueSum;
  if( this->VerifyNumberOfValidPoints( initialValueSum, derivative ) )
    {
    MeasureType valueSum( initialValueSum );
    //Threshold based on percentile
    size_t last_index = this->GetNumberOfValidPoints();
    if( m_Percentile < 100)
      {
      std::sort( values.begin(), values.end(), [](ValueType a, ValueType b)
        { return a.first < b.first ? true : false; });

      last_index = (this->GetNumberOfValidPoints() * this->m_Percentile) / 100;
      }
    if( last_index < 1 )
      {
      itkExceptionMacro( "Percentile too small. No valid point in selected percentile." );
      }

    //Accumulate percentile thresholded derivatives
    DerivativeType localTransformDerivative( this->GetNumberOfLocalParameters() );
    localTransformDerivative.Fill(0);
    const VirtualVectorContainer &virtualTransformedPointSet = 
      this->GetVirtualTransformedPointSet()->GetPoints()->CastToSTLConstContainer();  
    unsigned int nValidDistances = 0;
    for( typename std::vector<ValueType>::iterator it = values.begin();
         it < values.begin() + last_index; it++ )
      {
      ValueType &el = *it;
      PointIdentifier pointIndex = el.second.index;
      //Threshold based on distance
      if( el.first < this->m_DistanceCutoff )
        {
        valueSum += el.first;
        ++nValidDistances;

        if( this->HasLocalSupport() || this->GetCalculateValueAndDerivativeInTangentSpace() )
          {
          localTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );
          }
        for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
          {
          localTransformDerivative[par] += el.second.derivative[par];
          }

        if( this->HasLocalSupport() ||  this->GetCalculateValueAndDerivativeInTangentSpace() )
          {
          if( this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
            {
            this->StorePointDerivative( virtualTransformedPointSet[pointIndex], localTransformDerivative, derivative);
            }
          else
            {
            for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
              {
              derivative[this->GetNumberOfLocalParameters() * pointIndex + par] = localTransformDerivative[par];
              }
            }
          }
        }
      }
    if( nValidDistances < 1 )
      {
      itkExceptionMacro( "DistanceCutoff too small. No valid point in selected percentile." );
      }

    valueSum /= nValidDistances;
    if( !this->HasLocalSupport() && !this->GetCalculateValueAndDerivativeInTangentSpace() )
      {
      derivative = localTransformDerivative / nValidDistances;
      }
    calculatedValue = valueSum;
    this->m_Value = valueSum;
    }
  //clock2.Stop();

  //std::cout << "Collecting Derivatives Time: " << clock1.GetTotal() << std::endl;
  //std::cout << "Accumulating Derivatives Time: " << clock2.GetTotal() << std::endl;
}


/** PrintSelf */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Distance cuttoff: " << this->m_DistanceCutoff << std::endl;
  os << indent << "Percentile : " << this->m_Percentile << std::endl;
}
} // end namespace itk

#endif
