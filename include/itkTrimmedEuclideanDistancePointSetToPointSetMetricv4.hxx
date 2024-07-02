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
#ifndef itkTrimmedEuclideanDistancePointSetToPointSetMetricv4_hxx
#define itkTrimmedEuclideanDistancePointSetToPointSetMetricv4_hxx

#include "itkTrimmedEuclideanDistancePointSetToPointSetMetricv4.h"
#include "itkIdentityTransform.h"
#include "itkCompensatedSummation.h"

#include "itkMersenneTwisterRandomVariateGenerator.h"
#include <random>

namespace itk
{

/** Constructor */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::TrimmedEuclideanDistancePointSetToPointSetMetricv4()
{
  this->m_Percentile = 100;
  this->m_DistanceCutoff = NumericTraits<TInternalComputationValueType>::max();
  this->m_SamplingRate = 1.0;
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


  struct PointDerivativeStorage
    {
    MeasureType value;
    FixedPointIdentifier index;
    DerivativeType derivative;
    };

  using FixedTransformedVectorContainer = typename FixedPointsContainer::STLContainerType;
  using VirtualVectorContainer =  typename VirtualPointsContainer::STLContainerType;

  this->InitializeForIteration();

  derivative.SetSize( this->GetNumberOfParameters() );
  if( ! this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
    {
    derivative.SetSize( PointDimension * this->GetFixedTransformedPointSet()->GetNumberOfPoints() );
    }
  derivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  if( this->GetVirtualTransformedPointSet()->GetNumberOfPoints() !=
      this->GetFixedTransformedPointSet()->GetNumberOfPoints() )
    {
    itkExceptionMacro( "Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet." );
    }

  //Collect derviatives at each point
  std::vector< PointDerivativeStorage > values( this->GetFixedTransformedPointSet()->GetNumberOfPoints() );
  for(int i=0; i < values.size(); i++)
    {
    values[i].value = NumericTraits<MeasureType>::max();
    values[i].index = i;
    values[i].derivative.SetSize( this->GetNumberOfLocalParameters() );
    values[i].derivative.Fill(0);
    }
  {
  const VirtualVectorContainer &virtualTransformedPointSet =
    this->GetVirtualTransformedPointSet()->GetPoints()->CastToSTLConstContainer();
  const FixedTransformedVectorContainer &fixedTransformedPointSet =
    this->GetFixedTransformedPointSet()->GetPoints()->CastToSTLConstContainer();

  PointIdentifierRanges ranges = this->CreateRanges();
  std::function< void(int) > collectNeighborhoodValues =
       [ &values, this, &calculateValue, &derivative, &ranges,
         &virtualTransformedPointSet, &fixedTransformedPointSet]
       (int rangeIndex)
    {

    //The Jacobian might need to be alloacted per Thread for transforms with a large number of parmaters
    MovingTransformJacobianType jacobian( MovingPointDimension, this->GetNumberOfLocalParameters() );
    MovingTransformJacobianType jacobianCache;

    PixelType pixel;
    NumericTraits<PixelType>::SetLength( pixel, 1 );
    LocalDerivativeType pointDerivative;
    MeasureType pointValue = NumericTraits<MeasureType>::ZeroValue();

    //One generator per thread to avoid resource sharing
    using GeneratorType = itk::Statistics::MersenneTwisterRandomVariateGenerator;
    GeneratorType::Pointer generator = GeneratorType::New();
    generator->Initialize();

    /* Verify the virtual point is in the virtual domain.
     * If user hasn't defined a virtual space, and the active transform is not
     * a displacement field transform type, then this will always return true. */
    for(int index = ranges[rangeIndex].first; index<ranges[rangeIndex].second; index++)
      {

      if(this->m_SamplingRate < 1.0 )
        {
        if( generator->GetVariate() > this->m_SamplingRate )
          {
          continue;
          }
        }

      if( !this->IsInsideVirtualDomain( virtualTransformedPointSet[index] ) )
        {
        return;
        }
      if( this->m_UsePointSetData )
        {
        bool doesPointDataExist = this->GetFixedTransformedPointSet()->GetPointData( index, &pixel );
        if( ! doesPointDataExist )
          {
          itkExceptionMacro( "The corresponding data for point " << index << ") does not exist." );
          }
        }
      if( calculateValue || m_Percentile < 100 ||
          m_DistanceCutoff < NumericTraits<TInternalComputationValueType>::max() )
        {
        this->GetLocalNeighborhoodValueAndDerivative(
                fixedTransformedPointSet[index], pointValue, pointDerivative, pixel );
        values[index].value = pointValue;
        }
      else
        {
        pointDerivative = this->GetLocalNeighborhoodDerivative( fixedTransformedPointSet[index], pixel );
        }

      if( !this->GetCalculateValueAndDerivativeInTangentSpace() )
        {
        this->GetMovingTransform()->ComputeJacobianWithRespectToParametersCachedTemporaries(
                                       virtualTransformedPointSet[index], jacobian, jacobianCache );
        for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
          {
          for( DimensionType d = 0; d < PointDimension; ++d )
            {
            values[index].derivative[par] += jacobian(d, par) * pointDerivative[d];
            }
          }
        }
      else
        {
        for( DimensionType d = 0; d < PointDimension; ++d )
          {
          values[index].derivative[d] = pointDerivative[d];
          }
        }
      }
    };

  MultiThreaderBase::New()->ParallelizeArray( (FixedPointIdentifier) 0,
                          (FixedPointIdentifier) ranges.size(),
                           collectNeighborhoodValues, nullptr );
  }//End collection of derivatives

  // Threshold and sum derivaties and values
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
      std::sort( values.begin(), values.end(), [](PointDerivativeStorage a, PointDerivativeStorage b)
        { return a.value < b.value ? true : false; });

      last_index = ( values.size() * this->m_Percentile * this->m_SamplingRate ) / 100;
      }
    if( last_index < 1 )
      {
      itkExceptionMacro( "Percentile too small. No valid point in selected percentile." );
      }

    //Accumulate percentile thresholded derivatives and threshold by distance
    DerivativeType localTransformDerivative( this->GetNumberOfLocalParameters() );
    localTransformDerivative.Fill(0);
    const VirtualVectorContainer &virtualTransformedPointSet =
      this->GetVirtualTransformedPointSet()->GetPoints()->CastToSTLConstContainer();
    unsigned int nValidDistances = 0;
    for( int valueIndex=0; valueIndex < last_index; valueIndex++)
      {
      PointDerivativeStorage &el = values[valueIndex];
      PointIdentifier pointIndex = el.index;
      //Threshold based on distance
      if( el.value < this->m_DistanceCutoff )
        {
        valueSum += el.value;
        ++nValidDistances;

        if( this->HasLocalSupport() || this->GetCalculateValueAndDerivativeInTangentSpace() )
          {
          localTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );
          }
        for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
          {
          localTransformDerivative[par] += el.derivative[par];
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
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
const
typename TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PointIdentifierRanges
TrimmedEuclideanDistancePointSetToPointSetMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::CreateRanges() const
{
  PointIdentifier nPoints = this->m_FixedTransformedPointSet->GetNumberOfPoints();
  PointIdentifier nWorkUnits = MultiThreaderBase::New()->GetNumberOfWorkUnits();
  if( nWorkUnits > nPoints )
    {
    nWorkUnits = 1;
    }
  PointIdentifier startRange = 0;
  PointIdentifierRanges ranges;
  for(PointIdentifier p=1; p < nWorkUnits; ++p)
    {
    PointIdentifier endRange = (p * nPoints) / (double) nWorkUnits;
    ranges.push_back( PointIdentifierPair(startRange, endRange) );
    startRange = endRange;
    }
  ranges.push_back( PointIdentifierPair(startRange, nPoints) );
  return ranges;
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
