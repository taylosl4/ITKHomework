cmake_minimum_required(VERSION 2.8)

project (GeneralTemplate)

find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

add_executable (GeneralTemplate GeneralTemplate.cxx)

target_link_libraries (GeneralTemplate ${ITK_LIBRARIES})
