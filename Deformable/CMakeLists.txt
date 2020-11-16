cmake_minimum_required(VERSION 2.8)

project (DeformRegistration)

find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

add_executable (DeformRegistration DeformRegistration.cxx)

target_link_libraries (DeformRegistration ${ITK_LIBRARIES})
