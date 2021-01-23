include_directories(.)

add_library(Seldon SHARED Seldon.cxx ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})