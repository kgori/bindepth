cmake_minimum_required(VERSION 3.15)
project(BinDepth VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(ExternalProject)

find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBDEFLATE REQUIRED IMPORTED_TARGET libdeflate)

set(HTSLIB_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/htslib")
set(HTSLIB_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/htslib")
set(HTSLIB_INCLUDE_DIR "${HTSLIB_INSTALL_DIR}/include")
set(HTSLIB_STATIC_LIB "${HTSLIB_INSTALL_DIR}/lib/libhts.a")

ExternalProject_Add(
  htslib
  SOURCE_DIR ${HTSLIB_SOURCE_DIR}
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cd ${HTSLIB_SOURCE_DIR} && autoreconf -i
  COMMAND ${HTSLIB_SOURCE_DIR}/configure
          --prefix=${HTSLIB_INSTALL_DIR}
          --disable-shared
          --disable-s3
          --disable-gcs
          --disable-libcurl
          --disable-plugins
  BUILD_COMMAND $(MAKE)
  INSTALL_COMMAND $(MAKE) install
  BUILD_BYPRODUCTS ${HTSLIB_STATIC_LIB}
)

add_library(htslib::htslib STATIC IMPORTED GLOBAL)
file(MAKE_DIRECTORY ${HTSLIB_INCLUDE_DIR})
set_target_properties(htslib::htslib PROPERTIES
  IMPORTED_LOCATION ${HTSLIB_STATIC_LIB}
  INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${HTSLIB_INCLUDE_DIR}>"
)

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)

add_dependencies(htslib::htslib htslib)


add_executable(BinDepth main.cpp)
target_include_directories(BinDepth SYSTEM PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/CLI11")
target_link_libraries(BinDepth PRIVATE htslib::htslib ZLIB::ZLIB BZip2::BZip2 LibLZMA::LibLZMA)
target_compile_definitions(BinDepth PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:DEBUG_BUILD>)
if(UNIX AND NOT APPLE)
	target_link_libraries(BinDepth PRIVATE ${LIBDEFLATE_LIBRARIES})
	target_link_options(BinDepth PRIVATE -static-libstdc++ -static-libgcc)
endif()
