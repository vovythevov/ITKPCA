
#######################
#######################
#
#  Define Executables
#
#######################
#######################

FIND_PACKAGE(ITK REQUIRED COMPONENTS PrincipalComponentsAnalysis)

include_directories( ${PrincipalComponentsAnalysis_SOURCE_DIR}/Source )

ADD_EXECUTABLE(VectorKernelPCA VectorKernelPCA.cxx )
TARGET_LINK_LIBRARIES(VectorKernelPCA ${ITK_LIBRARIES})

