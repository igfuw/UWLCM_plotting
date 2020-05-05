cmake_minimum_required(VERSION 2.8.6) # not sure, simply copied from UWLCM
cmake_policy(SET CMP0028 NEW)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS iostreams program_options filesystem REQUIRED)

############################################################################################
#libcloudphxx
#find_package(libcloudph++)
#if (NOT libcloudph++_FOUND)
#  message(FATAL_ERROR "libclouphxx++ not found")
#endif()

############################################################################################
# UWLCM
find_package(UWLCM)
if (NOT UWLCM_FOUND)
  message(FATAL_ERROR "UWLCM not found")
endif()

if (NOT TARGET uwlcm::uwlcm)
  message(FATAL_ERROR "uwlcm::uwlcm target not found")
endif()

############################################################################################
# Blitz++
find_package(PkgConfig)
pkg_check_modules(BLITZ QUIET blitz>=0.10)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "Blitz++ library not found or Blitz++ version requirement not met (>=0.10)

* To insall Blitz++, please try:
*   Debian/Ubuntu: sudo apt-get install libblitz0-dev
*   Fedora: sudo yum install blitz-devel
*   Homebrew: brew install blitz
  ")
endif()
############################################################################################


foreach(exe drawbicyc plot_series_compare average)
  add_executable(${exe} ${exe}.cpp)
  # pointing to drawbicyc headers
  target_include_directories(${exe}
    PUBLIC
      $<INSTALL_INTERFACE:include>
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  )
  # enabling c++14
  #target_compile_features(${exe} PUBLIC cxx_std_14)
  #enabling c++11
  target_compile_features(${exe} PUBLIC cxx_std_11)
  target_link_libraries(${exe} PRIVATE uwlcm::uwlcm_includes)
  target_link_libraries(${exe} PRIVATE ${HDF5_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${exe} PRIVATE ${BLITZ_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${BLITZ_INCLUDE_DIRS})
  if(TARGET Boost::iostreams AND TARGET Boost::program_options AND TARGET Boost::filesystem)
    message("linkiing boost targets")
    target_link_libraries(${exe} PRIVATE Boost::iostreams Boost::program_options Boost::filesystem) 
  else()
    message("linkiing boost without targets")
    # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
    target_link_libraries(${exe} PRIVATE ${Boost_LIBRARIES})
    target_include_directories(${exe} PUBLIC ${Boost_INCLUDE_DIRS})
  endif()
endforeach()


############################################################################################
# installation

install(
  DIRECTORY
    include/drawbicyc
  DESTINATION
    include
)

install(TARGETS drawbicyc EXPORT drawbicyc-targets DESTINATION bin)
install(EXPORT drawbicyc-targets DESTINATION share/drawbicyc)

# generate and install a config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
  drawbicyc-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/drawbicyc-config.cmake"
  INSTALL_DESTINATION
    share/drawbicyc
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/drawbicyc-config.cmake"
  DESTINATION
    share/drawbicyc
)
############################################################################################



