# CMakeList.txt: CMake-Projekt für "ColDAEpp". Schließen Sie die Quelle ein, und definieren Sie
# projektspezifische Logik hier.
#
cmake_minimum_required (VERSION 3.8)

project ("ColDAEpp")

add_executable (ColDAEpp "ColDAEpp.cpp" "ColDAEpp.hpp" "linpack_d.cpp" "linpack_d.hpp" "blas0.cpp" "blas0.hpp" "blas1_d.cpp" "blas1_d.hpp" )

# TODO: Fügen Sie bei Bedarf Tests hinzu, und installieren Sie Ziele.
