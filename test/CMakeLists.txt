itk_module_test()

set(TrimmmedPointSetRegistrationTests
  itkMinimalStandardRandomVariateGeneratorTest.cxx
  itkMyFilterTest.cxx
  )

CreateTestDriver(TrimmmedPointSetRegistration "${TrimmmedPointSetRegistration-Test_LIBRARIES}" "${TrimmmedPointSetRegistrationTests}")

itk_add_test(NAME itkMinimalStandardRandomVariateGeneratorTest
  COMMAND TrimmmedPointSetRegistrationTestDriver itkMinimalStandardRandomVariateGeneratorTest
  )

itk_add_test(NAME itkMyFilterTest
  COMMAND TrimmmedPointSetRegistrationTestDriver
    --compare
    ${ITK_TEST_OUTPUT_DIR}/itkMyFilterTestOutput.mha
    DATA{Baseline/itkMyFilterTestOutput.mha}
  itkMyFilterTest
    ${ITK_TEST_OUTPUT_DIR}/itkMyFilterTestOutput.mha
  )
