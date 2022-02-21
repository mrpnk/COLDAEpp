# CMakeList.txt: CMake-Projekt für "ColDAEpp". Schließen Sie die Quelle ein, und definieren Sie
# projektspezifische Logik hier.
#
cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

project ("ColDAEpp")

#target_compile_features(ColDAEpp PRIVATE cxx_std_17)

file(GLOB_RECURSE SOURCES src/*.cpp)

add_executable (ColDAEpp ${SOURCES})


# TODO: Fügen Sie bei Bedarf Tests hinzu, und installieren Sie Ziele.
