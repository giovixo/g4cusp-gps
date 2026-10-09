#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(cusp-activation)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Sources: one directory per part, in both src/ and include/
#   common/          geometry and physics list, shared by the two programs
#   activation/      user actions of cusp-activation (nuclide production)
#   postactivation/  user actions of cusp-postactivation (decay of the nuclides)
# NB: headers are listed so they show up in IDEs.
# CONFIGURE_DEPENDS re-runs the globs at build time, so added/removed files are picked up
#
file(GLOB common_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/common/*.cc)
file(GLOB common_headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/common/*.hh)
file(GLOB activation_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/activation/*.cc)
file(GLOB activation_headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/activation/*.hh)
file(GLOB postactivation_sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/postactivation/*.cc)
file(GLOB postactivation_headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/postactivation/*.hh)

#----------------------------------------------------------------------------
# Shared code, built once as a static library
#
add_library(cusp-common STATIC ${common_sources} ${common_headers})
target_include_directories(cusp-common PUBLIC ${PROJECT_SOURCE_DIR}/include/common)
target_link_libraries(cusp-common PUBLIC ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Add the executables, and link them to the shared code and the Geant4 libraries
#
add_executable(cusp-activation cusp-activation.cc ${activation_sources} ${activation_headers})
target_include_directories(cusp-activation PRIVATE ${PROJECT_SOURCE_DIR}/include/activation)
target_link_libraries(cusp-activation cusp-common)

add_executable(cusp-postactivation cusp-postactivation.cc ${postactivation_sources} ${postactivation_headers})
target_include_directories(cusp-postactivation PRIVATE ${PROJECT_SOURCE_DIR}/include/postactivation)
target_link_libraries(cusp-postactivation cusp-common)

#----------------------------------------------------------------------------
# Copy the macro and GDML files to the build directory at every build
# (the program reads them from the working directory)
#
file(GLOB RUNTIME_FILES CONFIGURE_DEPENDS
     ${PROJECT_SOURCE_DIR}/macros/*.mac
     ${PROJECT_SOURCE_DIR}/gdml-mass-model/*.xml
     ${PROJECT_SOURCE_DIR}/gdml-mass-model/*.gdml)
add_custom_target(copy_runtime_files ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${RUNTIME_FILES} ${PROJECT_BINARY_DIR}
                  DEPENDS ${RUNTIME_FILES}
                  COMMENT "Copying macro and GDML files to ${PROJECT_BINARY_DIR}")
add_dependencies(cusp-activation copy_runtime_files)
add_dependencies(cusp-postactivation copy_runtime_files)

#----------------------------------------------------------------------------
# Install the executables to 'bin/Darwin-g++' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS cusp-activation cusp-postactivation DESTINATION bin/Darwin-g++)
