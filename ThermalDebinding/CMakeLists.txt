##
#  CMake script for the ThermalDebinding program:
##

# Set the name of the project and target:
SET(TARGET "ThermalDebinding")

SET(TARGET_SRC
  ${TARGET}.cc Parameters.cc Material.cc Time.cc
  )

# Define the output that should be cleaned:
SET(CLEAN_UP_FILES *.vtu)

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.3.0)

FIND_PACKAGE(deal.II 9.4.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

IF(NOT DEAL_II_WITH_SUNDIALS) # keep in one line
  MESSAGE(FATAL_ERROR "
Error! This program requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_SUNDIALS = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options:
    DEAL_II_WITH_SUNDIALS = ${DEAL_II_WITH_SUNDIALS}
This conflicts with the requirements."
    )
ENDIF()

add_compile_options(-g)

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
