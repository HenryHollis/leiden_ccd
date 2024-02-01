cmake_minimum_required(VERSION 3.27)
project(leiden_ccd)

set(CMAKE_CXX_STANDARD 17)

find_package(igraph REQUIRED)

add_library(libleidenalg
        src/GraphHelper.cpp
        src/ModularityVertexPartition.cpp
        src/MutableVertexPartition.cpp
        src/Optimiser.cpp
        include/GraphHelper.h
        include/MutableVertexPartition.h
        include/ModularityVertexPartition.h
        include/Optimiser.h
        src/ccdModularityVertexPartition.cpp
        include/ccdModularityVertexPartition.h
        src/ccd_utils.cpp
        include/ccd_utils.h
)


# Generate export header include file
include(GenerateExportHeader)
generate_export_header(libleidenalg
        STATIC_DEFINE LEIDENALG_STATIC
        EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/libleidenalg_export.h
)

# Make sure we install the export header
target_sources(libleidenalg
        PUBLIC
        FILE_SET headers
        TYPE HEADERS
        BASE_DIRS
        ${PROJECT_BINARY_DIR}/include/
        FILES
        ${PROJECT_BINARY_DIR}/include/libleidenalg_export.h
)

# Make sure we can find include files
target_include_directories(
        libleidenalg
        PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)


target_link_libraries (libleidenalg PUBLIC igraph::igraph)

#For adding example.cpp with Main()
add_executable(example tests/example.cpp)
add_executable(example2 tests/example2.cpp)
add_executable(test_calcCorMat tests/test_calcCorMat.cpp)
add_executable(test_calcCCDsimple tests/test_calcCCDsimple.cpp)
add_executable(simple_igraph_creation tests/simple_igraph_creation.cpp)
add_executable(Multiplex_optimization_example tests/Multiplex_optimization_example.cpp)

target_link_libraries(example libleidenalg)
target_link_libraries(example2 libleidenalg)
target_link_libraries(test_calcCorMat libleidenalg)
target_link_libraries(test_calcCCDsimple libleidenalg)
target_link_libraries(simple_igraph_creation libleidenalg)
target_link_libraries(Multiplex_optimization_example libleidenalg)
