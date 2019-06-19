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

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})


#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PBNRR)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
