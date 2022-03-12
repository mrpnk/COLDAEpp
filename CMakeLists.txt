# CMakeList.txt: CMake-Projekt für "ColDAEpp".
cmake_minimum_required (VERSION 3.8)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:fast")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")

project ("ColDAEpp")

# define the 'sources' to be all the .cpp files in this directory and all subdirectories
file(GLOB_RECURSE SOURCES test/linpack/*.cpp)


# Here come the tests

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


enable_testing()

add_executable (testAll test/test.cpp ${SOURCES})
target_link_libraries(testAll gtest_main)


# set some options for the project
#target_compile_options(ColDAEpp PRIVATE /W4)
target_compile_features(testAll PRIVATE cxx_std_20)



include(GoogleTest)
gtest_discover_tests(testAll)
