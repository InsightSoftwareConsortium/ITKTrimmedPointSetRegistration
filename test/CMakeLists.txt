itk_module_test()

set(TrimmedPointSetRegistrationTests
        itkTrimmedEuclideanDistancePointSetRegistrationTest.cxx
   )

CreateTestDriver(TrimmedPointSetRegistration "${TrimmedPointSetRegistration-Test_LIBRARIES}" "${TrimmedPointSetRegistrationTests}")

itk_add_test(NAME itkTrimmedEuclideanDistancePointSetRegistrationTest
        COMMAND TrimmedPointSetRegistrationTestDriver itkTrimmedEuclideanDistancePointSetRegistrationTest
  )
