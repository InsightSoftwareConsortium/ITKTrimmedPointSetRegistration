itk_module_test()

set(TrimmedPointSetRegistrationTests
      itkTrimmedPointSetRegistrationTest.cxx
   )

CreateTestDriver(TrimmedPointSetRegistration "${TrimmedPointSetRegistration-Test_LIBRARIES}" "${TrimmedPointSetRegistrationTests}")

itk_add_test(NAME itkPointSetRegistrationTest
  COMMAND TrimmedPointSetRegistrationTestDriver itkTrimmedPointSetRegistrationTest
  )
