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
  this->InitializeForIteration();

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  PointIdentifier numberOfFixedPoints = this->GetFixedTransformedPointSet()->GetNumberOfPoints();
  if( this->GetVirtualTransformedPointSet()->GetNumberOfPoints() != numberOfFixedPoints )
    {
    itkExceptionMacro("Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet.");
    }

  std::vector<MeasureType> distances(numberOfFixedPoints, NumericTraits<MeasureType>::max());

  std::function< void(PointsConstIterator) > collectNeighborhoodValues =
          [&distances, this](PointsConstIterator It)
    {
    FixedPointIdentifier index = It.Index();
    FixedPointType virtualTransformedPoint = this->GetVirtualTransformedPointSet()->GetPoint( index);
    FixedPointType fixedTransformedPoint =  this->GetFixedTransformedPointSet()->GetPoint( index);
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

  std::for_each( this->GetFixedTransformedPointSet()->GetPoints()->Begin(),
                 this->GetFixedTransformedPointSet()->GetPoints()->End(), collectNeighborhoodValues );

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
  using ValueType = std::pair<MeasureType, DerivativeType>;

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

 // std::function< void(PointsConstIterator) > collectNeighborhoodValues =
 //         [&values, this, &calculateValue, &derivative](PointsConstIterator It)
 PointsConstIterator virtualIt = this->GetVirtualTransformedPointSet()->GetPoints()->Begin();
 PointsConstIterator It = this->GetFixedTransformedPointSet()->GetPoints()->Begin();
 PointsConstIterator end = this->GetFixedTransformedPointSet()->GetPoints()->End();
 for( ;It != end; ++It, ++virtualIt )
   {

   thread_local MovingTransformJacobianType jacobian( MovingPointDimension, this->GetNumberOfLocalParameters() );
   thread_local MovingTransformJacobianType jacobianCache;
   thread_local DerivativeType localTransformDerivative( this->GetNumberOfLocalParameters() );

   localTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

   MeasureType pointValue = NumericTraits<MeasureType>::ZeroValue();
   LocalDerivativeType pointDerivative;

   FixedPointIdentifier index = It.Index();
   const FixedPointType &virtualTransformedPoint = virtualIt.Value();
   const FixedPointType &fixedTransformedPoint =  It.Value();

   /* Verify the virtual point is in the virtual domain.
    * If user hasn't defined a virtual space, and the active transform is not
    * a displacement field transform type, then this will always return true. */
   values[index].first = 0;
   if( this->IsInsideVirtualDomain( virtualTransformedPoint ) )
     {
     PixelType pixel;
     NumericTraits<PixelType>::SetLength( pixel, 1 );
     if( this->m_UsePointSetData )
       {
       bool doesPointDataExist = this->GetFixedTransformedPointSet()->GetPointData( It.Index(), &pixel );
       if( ! doesPointDataExist )
         {
         itkExceptionMacro( "The corresponding data for point " << It.Value() << " (pointId = " << It.Index() << ") does not exist." );
         }
       }
     if( calculateValue )
       {
       this->GetLocalNeighborhoodValueAndDerivative( fixedTransformedPoint, pointValue, pointDerivative, pixel );
       values[index].first = pointValue;
       }
     else
       {
       pointDerivative = this->GetLocalNeighborhoodDerivative( fixedTransformedPoint, pixel );
       }

     if( !this->GetCalculateValueAndDerivativeInTangentSpace() )
       {
       this->GetMovingTransform()->
          ComputeJacobianWithRespectToParametersCachedTemporaries( virtualTransformedPoint,
                                                                   jacobian,
                                                                   jacobianCache );
       for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
         {
         for( DimensionType d = 0; d < PointDimension; ++d )
           {
           localTransformDerivative[par] += jacobian(d, par) * pointDerivative[d];
           }
         }
       }
     else
       {
       for( DimensionType d = 0; d < PointDimension; ++d )
         {
         localTransformDerivative[d] += pointDerivative[d];
         }
       }
      // For local-support transforms, store the per-point result
     if( this->HasLocalSupport() || this->GetCalculateValueAndDerivativeInTangentSpace() )
       {
       if( this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
         {
         this->StorePointDerivative( virtualTransformedPoint, localTransformDerivative, derivative );
         }
       else
         {
         for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
           {
           derivative[this->GetNumberOfLocalParameters() * It.Index() + par] = localTransformDerivative[par];
           }
         }
       }
     else
       {
       values[index].second = localTransformDerivative;
       }
     }
   };

  //std::for_each( this->GetFixedTransformedPointSet()->GetPoints()->Begin(),
  //               this->GetFixedTransformedPointSet()->GetPoints()->End(), collectNeighborhoodValues );


  // `valueSum` default value is set to max in `VerifyNumberOfValidPoints`
  // if there is no valid point.
  MeasureType valueSum = 0.0;
  if( this->VerifyNumberOfValidPoints( valueSum, derivative ) )
    {
    size_t last_index = this->GetNumberOfValidPoints();
    if( m_Percentile < 100)
      {
      std::sort( values.begin(), values.end(), [](ValueType a, ValueType b)
        { return a.first < b.first ? true : false; });

      last_index = this->GetNumberOfValidPoints() * m_Percentile / 100;
      }
    if( last_index < 1 )
      {
      itkExceptionMacro( "Percentile too small. No valid point in selected percentile." );
      }

    unsigned int nDistances = 0;
    for( typename std::vector<ValueType>::iterator it = values.begin();
         it < values.begin() + last_index; it++ )
      {
      ValueType &el = *it;
      if( el.first < m_DistanceCutoff )
        {
        valueSum += el.first;
        ++nDistances;
        if( ! this->HasLocalSupport() && ! this->GetCalculateValueAndDerivativeInTangentSpace() )
          {
          for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
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
  this->m_Value = valueSum;
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
