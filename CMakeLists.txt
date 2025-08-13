include(CPack)

cmake_minimum_required(VERSION 3.10.0)
project(SimplePKPDProject VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(SimplePKPDProject 
    main.cpp
    pk_models/one_compartment.cpp
)

target_compile_features(SimplePKPDProject PRIVATE cxx_std_20)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
