cmake_minimum_required(VERSION 3.5)

project(PBNRR)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/aangelos28/PBNRR")
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Angelos Angelopoulos (CRTC), Fotis Drakopoulos (CRTC), Yixun Liu (CRTC), Andriy Kot (CRTC), Andrey Fedorov (SPL B&W Harvard), Olivier Clatz (Asclepios INRIA), Nikos Chrisochoides (CRTC)")
set(EXTENSION_DESCRIPTION "This Slicer extension non-rigid registers a moving to a fixed MRI using a linear homogeneous bio-mechanical model to compute a dense deformation field that defines a transformation for every point in the fixed to the moving image.")
set(EXTENSION_ICONURL "https://raw.github.com/aangelos28/PBNRR/master/PBNRR.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/PBNRR/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any
set(EXTENSION_BUILD_SUBDIRECTORY inner-build)

set(SUPERBUILD_TOPLEVEL_PROJECT inner)

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
mark_as_superbuild(Slicer_DIR)

find_package(Git REQUIRED)
mark_as_superbuild(GIT_EXECUTABLE)

#-----------------------------------------------------------------------------
# SuperBuild setup
option(${EXTENSION_NAME}_SUPERBUILD "Build ${EXTENSION_NAME} and the projects it depends on." ON)
mark_as_advanced(${EXTENSION_NAME}_SUPERBUILD)
if(${EXTENSION_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PBNRR)
## NEXT_MODULE

#-----------------------------------------------------------------------------
set(EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS)
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${ITKFEM_DIR};ITKFEM;RuntimeLibraries;/")
set(${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS "${EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS}" CACHE STRING "List of external projects to install" FORCE)

#-----------------------------------------------------------------------------
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS}")
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
