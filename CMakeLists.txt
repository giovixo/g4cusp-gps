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
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
# CONFIGURE_DEPENDS re-runs the glob at build time, so added/removed files are picked up
#
file(GLOB sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(cusp-activation cusp-activation.cc ${sources} ${headers})
target_link_libraries(cusp-activation ${Geant4_LIBRARIES})

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

#----------------------------------------------------------------------------
# Install the executable to 'bin/Darwin-g++' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS cusp-activation DESTINATION bin/Darwin-g++)
