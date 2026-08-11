cmake_minimum_required(VERSION 3.15)

option(STANDALONE_BUILD "Build for system package manager" OFF)
set(STANDALONE_BUILD ${STANDALONE_BUILD} CACHE BOOL "Pass standalone build flag" FORCE)

# Set CMP0167 to NEW to use the modern FindBoost behavior
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
project(pyjeo LANGUAGES CXX)

# Must be at the start of the file
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 1. RPATH setup (Crucial for the "Build" case, ignored for the "System" case)
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif()
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

find_package(Git REQUIRED)

set(CMAKE_DISABLE_FIND_PACKAGE_Boost ON)
set(Boost_NO_BOOST_CMAKE ON CACHE BOOL "Disable Boost Config mode" FORCE)

# Voeg dit toe zodat CMake de libraries direct op de juiste plek zoekt:
set(BOOST_ROOT "/usr" CACHE PATH "Set Boost root" FORCE)
set(BOOST_INCLUDEDIR "/usr/include" CACHE PATH "Set Boost include" FORCE)
set(BOOST_LIBRARYDIR "/usr/lib64" CACHE PATH "Set Boost libs" FORCE)

set(Python_EXECUTABLE "${Python_EXECUTABLE}" CACHE FILEPATH "Python executable" FORCE)
set(Python_INCLUDE_DIRS "${Python_INCLUDE_DIRS}" CACHE PATH "Python include dirs" FORCE)
set(Python_NumPy_INCLUDE_DIRS "${Python_NumPy_INCLUDE_DIRS}" CACHE PATH "NumPy include" FORCE)
set(Python_Development_Module_FOUND TRUE CACHE BOOL "Force Python dev module found" FORCE)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/miallib/CMakeLists.txt")
    message(STATUS "Submodule miallib missing. Retrieving miallib from git repository...")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --remote -- external/miallib
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/jiplib/CMakeLists.txt")
    message(STATUS "Submodule jiplib missing. Retrieving jiplib from git repository...")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --remote -- external/jiplib
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

message(STATUS "Looking for miallib")
if(STANDALONE_BUILD)
    find_package(miallib 1.2.1 REQUIRED) # Debian
else()
    find_package(miallib 1.2.1 QUIET NO_DEFAULT_PATH)
endif()

if(NOT miallib_FOUND)
    message(STATUS "miallib not found, building from source")
    add_subdirectory(external/miallib)
else()
    message(STATUS "Found system miallib")
endif()

message(STATUS "Looking for jiplib")
if(STANDALONE_BUILD)
    find_package(jiplib 1.2.1 REQUIRED) # Debian
else()
    find_package(jiplib 1.2.1 QUIET NO_DEFAULT_PATH)
endif()
if(NOT jiplib_FOUND)
    message(STATUS "jiplib not found, building from source")
    add_subdirectory(external/jiplib)
else()
    message(STATUS "Found system jiplib")
endif()

#find_library(JSONCPP_LIB_PATH NAMES jsoncpp PATHS ${JSONCPP_LIBRARY_DIRS})

if(NOT jiplib_FOUND)
  # 4. Define the Python Extension (SWIG)
  find_package(SWIG REQUIRED)
  find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
  include(${SWIG_USE_FILE})
  
  # Ensure we can find headers regardless of source
  target_include_directories(jiplib_python PRIVATE
      Python::Module
      #${PYTHON_INCLUDE_DIRS}
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/jiplib/src>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/miallib/src>
      #external/jiplib/src
      #external/miallib/src
  )
  
  # 5. Smart Linking
  # Link against the targets (CMake handles whether they are system or local)
  target_link_libraries(jiplib_python PRIVATE jiplib miallib Python::Module)
  #Ensure jiplib is compiled with a standard high enough for JsonCpp
  target_compile_features(jiplib PRIVATE cxx_std_11)
endif()
  
# 6. Installation
# If we built them locally, we install them into the wheel root
if(NOT miallib_FOUND)
    install(TARGETS miallib DESTINATION .)
endif()
if(NOT jiplib_FOUND)
    install(TARGETS jiplib DESTINATION .)
    install(TARGETS jiplib_python DESTINATION .)
    install(FILES $<TARGET_PROPERTY:jiplib_python,SWIG_SUPPORT_FILES_DIRECTORY>/jiplib.py DESTINATION .)
endif()

install(DIRECTORY pyjeo/ DESTINATION pyjeo)
