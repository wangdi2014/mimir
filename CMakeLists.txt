cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")


project(mimir)

file(GLOB SelectionSources "source/selection/*.h" "source/selection/*.cpp")
file(GLOB AppSources "source/app/*.cpp")


include_directories(${CMAKE_SOURCE_DIR}/eigen)
include_directories(${CMAKE_SOURCE_DIR}/source/selection)
include_directories(${CMAKE_SOURCE_DIR}/source/app)

add_library(selection STATIC ${SelectionSources})
add_executable(app source/app/app.cpp)

target_link_libraries(app selection)