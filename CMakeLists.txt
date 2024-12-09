cmake_minimum_required(VERSION 3.29)
project(FluidSimulation)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify compilers
set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
# Add Raylib
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # Don't build the supplied examples
set(BUILD_GAMES OFF CACHE BOOL "" FORCE)    # Don't build the supplied example games

FetchContent_Declare(
        raylib
        GIT_REPOSITORY "https://github.com/raysan5/raylib.git"
        GIT_TAG "master"
        GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(raylib)

# Debug Raylib
message(STATUS "Raylib source directory: ${raylib_SOURCE_DIR}")
message(STATUS "Raylib binary directory: ${raylib_BINARY_DIR}")

# Add sources
file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/sources/*.cpp")
set(PROJECT_INCLUDE "${CMAKE_CURRENT_LIST_DIR}/sources/")

# Find OpenMP
find_package(OpenMP REQUIRED COMPONENTS CXX)

if (OpenMP_FOUND)
    message(STATUS "Found OpenMP for C++: ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "Could not find OpenMP for C++")
endif()

# Declare executable
add_executable(${PROJECT_NAME} main.cpp
        SimulationParticle2D.cpp
        SimulationParticle2D.h
)

# Link OpenMP and raylib
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX raylib)

# Include project sources
target_sources(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_INCLUDE})

# Set asset path
target_compile_definitions(${PROJECT_NAME} PUBLIC ASSETS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/assets/")

# macOS-specific framework linking
if (APPLE)
    target_link_libraries(${PROJECT_NAME} PRIVATE
            "/System/Library/Frameworks/OpenGL.framework"
            "/System/Library/Frameworks/Cocoa.framework"
            "/System/Library/Frameworks/IOKit.framework"
            "/System/Library/Frameworks/CoreVideo.framework"
    )
endif()