cmake_minimum_required(VERSION 3.1...3.14)

if(${CMAKE_VERSION} VERSION_LESS 3.14)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.14)
endif()

# Force out-of-source build
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

## Set correct project version based on: https://semver.org/
# Given a version number MAJOR.MINOR.PATCH, increment the:
#
# MAJOR version when you make incompatible API changes,
# MINOR version when you add functionality in a backwards-compatible manner, and
# PATCH version when you make backwards-compatible bug fixes.

project(Unpacker2 VERSION 2.0.1
                  LANGUAGES CXX)
set(PROJECT_DESCRIPTION "Unpacker2 module")

message(STATUS "")
message(STATUS "    == ${PROJECT_NAME} Project configuration ==")
message(STATUS "")

set(CMAKE_CXX_STANDARD 11) # maybe not best solution to set it globally? upgrate minimal cmake version to 3.8 or find another way to set it only for target with interface
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Include cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Be nice and export compile commands by default, this is handy for clang-tidy
# and for other tools.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Helpful option enable build profiling to identify slowly compiling files
option(MEASURE_ALL "When enabled all commands will be passed through time command" OFF)
if(MEASURE_ALL)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "time")
endif()

################################################################################
## Install format hook to git
## Automatically installs hook to .git/hooks/pre-commit to format code with clang-format
option(INSTALL_HOOK "Install format hook to .git/hooks/pre-commit" ON)
if(INSTALL_HOOK)
  message(STATUS "Trying to install format hook..")
  set(INSTALL_HOOK ${CMAKE_CURRENT_SOURCE_DIR}/hooks/git-pre-commit-format install)
  execute_process(COMMAND ${INSTALL_HOOK})
endif()

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    Unpacker2ConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

#Configure Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58 REQUIRED) # boost component for header-only libraries

if(NOT TARGET Boost::boost)
    add_library(Boost::boost IMPORTED INTERFACE)
    set_property(TARGET Boost::boost PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::boost PROPERTY
        INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

# Configure ROOT
## Based on https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
find_package(ROOT CONFIG REQUIRED)
if(ROOT_FOUND)
  message(STATUS "ROOT (version ${ROOT_VERSION}) was found using ROOTConfig")
endif()

if(ROOT_VERSION VERSION_LESS "6.10")
  message(FATAL_ERROR "Only ROOT >= 6.10 is currently supported")
endif(ROOT_VERSION VERSION_LESS "6.10")

if(EXISTS "${ROOT_DIR}/modules/RootNewMacros.cmake")
  include("${ROOT_DIR}/modules/RootNewMacros.cmake")
  message(STATUS "Found RootNewMacros.cmake in: ${ROOT_DIR}/modules!")
elseif(EXISTS "${ROOT_DIR}/cmake/modules/RootNewMacros.cmake")
  include("${ROOT_DIR}/cmake/modules/RootNewMacros.cmake")
  message(STATUS "Found RootNewMacros.cmake in: ${ROOT_DIR}/cmake/modules!")
elseif(EXISTS "${ROOT_DIR}/../cmake/modules/RootNewMacros.cmake")
  include("${ROOT_DIR}/../cmake/modules/RootNewMacros.cmake")
  message(STATUS "Found RootNewMacros.cmake in: ${ROOT_DIR}/../cmake/modules!")
else()
  message("ERROR! Could not find RootNewMacros.cmake file! Aborting..")
  return(1)
endif()
# fix missing include directories property in versions < ROOT 6.12
set_target_properties(ROOT::Core PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")


if(NOT TARGET ROOT::Flags_CXX)
    add_library(ROOT::Flags_CXX IMPORTED INTERFACE)
    # Fix for ROOT_CXX_FLAGS not actually being a CMake list
    separate_arguments(ROOT_CXX_FLAGS)
    set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
        INTERFACE_COMPILE_OPTIONS ${ROOT_CXX_FLAGS})
endif(NOT TARGET ROOT::Flags_CXX)


# Add definitions
separate_arguments(ROOT_DEFINITIONS)
foreach(_flag ${ROOT_EXE_LINKER_FLAG_LIST})
    # Remove -D or /D if present
    string(REGEX REPLACE [=[^[-//]D]=] "" _flag ${_flag})
    set_property(TARGET ROOT::Flags APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${_flag})

    # This also fixes a bug in the linker flags
    string(REPLACE "-L " "-L" ROOT_EXE_LINKER_FLAGS "${ROOT_EXE_LINKER_FLAGS}")
    separate_arguments(ROOT_EXE_LINKER_FLAGS)
    set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES ${ROOT_EXE_LINKER_FLAGS})
endforeach()
set(DICTIONARY_NAME G__Unpacker2)
set(HEADERS_WITH_DICTIONARY_REQUIRED
 ${CMAKE_CURRENT_SOURCE_DIR}/include/EventIII.h
 ${CMAKE_CURRENT_SOURCE_DIR}/include/TDCChannel.h
  )

set(LINKDEF_NAME Unpacker2_LinkDef.h)

include_directories(ROOT_BUG) # Fix bug in NewMacros Root file
ROOT_GENERATE_DICTIONARY(${DICTIONARY_NAME} ${HEADERS_WITH_DICTIONARY_REQUIRED}
  OPTIONS -Iinclude
  LINKDEF ${LINKDEF_NAME})
## End of configuration of ROOT

# Set all sources
set(SOURCES src/Unpacker2.cc src/Unpacker2D.cc src/EventIII.cc src/TDCChannel.cc ${DICTIONARY_NAME}.cxx)

add_library(Unpacker2 SHARED ${SOURCES})
#TODO: maybe change namespace to JPetTomography?
add_library(Unpacker2::Unpacker2 ALIAS Unpacker2)
target_include_directories(Unpacker2
                           PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                           $<INSTALL_INTERFACE:include>)
target_link_libraries(Unpacker2 PUBLIC Boost::boost ROOT::Core ROOT::Hist ROOT::Tree ROOT::Flags_CXX)
set_target_properties(Unpacker2 PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(Unpacker2 PROPERTIES POSITION_INDEPENDENT_CODE ON) #-fPIC for shared library
set_target_properties(Unpacker2 PROPERTIES VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})

#include(CheckIPOSupported) # only avaiable in 3.9+
#check_ipo_supported(RESULT IPOresult)
#if(IPOresult)
#  set_target_properties(Unpacker2 PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
#endif()

include(GNUInstallDirs)
install(TARGETS Unpacker2
        EXPORT Unpacker2
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT Unpacker2
        FILE Unpacker2Targets.cmake
        NAMESPACE Unpacker2::
        DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
         )

configure_file(Unpacker2Config.cmake.in Unpacker2Config.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Unpacker2Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/Unpacker2ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
        )
        
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libUnpacker2_rdict.pcm"
              "${CMAKE_CURRENT_BINARY_DIR}/libUnpacker2.rootmap"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

        configure_file("${PROJECT_BINARY_DIR}/Unpacker2Config.cmake" "${PROJECT_BINARY_DIR}/cmake/Unpacker2Config.cmake" COPYONLY) 
        configure_file("${PROJECT_BINARY_DIR}/Unpacker2ConfigVersion.cmake" "${PROJECT_BINARY_DIR}/cmake/Unpacker2ConfigVersion.cmake" COPYONLY) 
export(EXPORT Unpacker2 FILE Unpacker2Targets.cmake NAMESPACE Unpacker2::)
add_custom_target(
  copy_target
  COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_BINARY_DIR}/Unpacker2Targets.cmake" "${PROJECT_BINARY_DIR}/cmake/Unpacker2Targets.cmake"
  DEPENDS "${PROJECT_BINARY_DIR}/Unpacker2Targets.cmake")
add_dependencies(Unpacker2 copy_target) 

#unit tests
option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Packaging support
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VENDOR "JPetTomography")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Package contains Unpacker2 module")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "grey")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
#TODO: add missing licence
#set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /*.dist*
    /.*build.*
    /\\\\.DS_Store
)

include(CPack)

################################################################################
# Wrap up of settings printed on build
message(STATUS "")
message(STATUS "    == Final overview for ${PROJECT_NAME} ==")
message(STATUS "Version:              ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
message(STATUS "Install prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Compiler:             ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}")
message(STATUS "  possible options: Debug Release RelWithDebInfo MinSizeRel")
message(STATUS "  set with ` cmake -DCMAKE_BUILD_TYPE=Debug .. `")
message(STATUS "")
