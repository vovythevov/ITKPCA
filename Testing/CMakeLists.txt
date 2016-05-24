SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
SET(TEST_DATA_ROOT ${CMAKE_SOURCE_DIR}/Data)
SET(TEMP ${CMAKE_BINARY_DIR}/Testing/Temporary)


#######################
#######################
#
#  Define Executables
#
#######################
#######################

ADD_EXECUTABLE(itkVectorKernelPCATest itkVectorKernelPCATest.cxx )
TARGET_LINK_LIBRARIES(itkVectorKernelPCATest ${ITK_LIBRARIES})

###################
###################
#
#  Define Tests
#
###################
###################

ADD_TEST(itkVectorKernelPCATest
  ${CXX_TEST_PATH}/itkVectorKernelPCATest
  ${TEST_DATA_ROOT}/PCATestSurface.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_01.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_02.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_03.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_04.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_05.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_06.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_07.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_08.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_09.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_10.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_11.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_12.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_13.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_14.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_15.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_16.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_17.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_18.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_19.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_20.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_21.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_22.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_23.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_24.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_25.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_26.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_27.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_28.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_29.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_30.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_31.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_32.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_33.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_34.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_35.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_36.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_37.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_38.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_39.vtk
  ${TEST_DATA_ROOT}/PCATestSurface_alpha0_40.vtk
  )

