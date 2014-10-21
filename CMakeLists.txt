
# Set the name of the project and target:
SET(TARGET "quantum_dot_problem")

SET(TARGET_SRC
  source/geometry_function.cc
  source/parameter_reader.cc
#
  source/dielectric_tensor.cc
  source/elastic_tensor.cc
  source/piezoelectric_tensor.cc
  source/polarelectric_tensor.cc
#
  source/base/group_symmetry.cc
  source/base/tensor_base.cc
#
  six_band_hole_model.cc
  piezoelectric_model.cc
  piezoelectric_coefficients.cc
#
  ${TARGET}.cc
  )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.2 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
