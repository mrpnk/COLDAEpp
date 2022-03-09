# CMakeList.txt: CMake-Projekt für "ColDAEpp".
cmake_minimum_required (VERSION 3.8)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:precise")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")

project ("ColDAEpp")

# define the 'sources' to be all the .cpp files in this directory and all subdirectories
file(GLOB_RECURSE SOURCES src/*.cpp)

# make the target 'ColDAEpp'
add_executable (ColDAEpp ${SOURCES})

# set some options for the project
target_compile_options(ColDAEpp PRIVATE /W4)
target_compile_features(ColDAEpp PRIVATE cxx_std_20)


# TODO: Fügen Sie bei Bedarf Tests hinzu, und installieren Sie Ziele.
