cmake_minimum_required(VERSION 2.8)

project (AffineRegistration)

find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

add_executable (AffineRegistration AffineRegistration.cxx)

target_link_libraries (AffineRegistration ${ITK_LIBRARIES})
