cmake_minimum_required(VERSION 3.10)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)
set(GCC_MIN_VERSION "9") # minimum GCC version that is tested in CI

project(AliStat)

set(CPACK_PACKAGE_VERSION "1.16")

###############################################################################
## file globbing ##############################################################
###############################################################################

file(GLOB_RECURSE sources      src/*.cpp src/*.h)

###############################################################################
## target definitions #########################################################
###############################################################################

add_executable(alistat ${sources})

# this lets me include files relative to the root source directory with a <> pair
target_include_directories(alistat PUBLIC src)

##################################################################
# Detect target platforms
##################################################################
if (WIN32)
    message("Target OS     : Windows")
    # build as static binary to run on most machines
    if (ALISTAT_FLAGS MATCHES "static")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif (APPLE)
    message("Target OS     : Mac OS X")
        if(OSX_NATIVE_ARCHITECTURE STREQUAL "arm64")
                add_definitions("--target=arm64-apple-macos10.5")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --target=arm64-apple-macos12.0.1")
        else()
                # to be compatible back to Mac OS X 10.7
        if (ALISTAT_FLAGS MATCHES "oldmac")
                add_definitions("-mmacosx-version-min=10.5")
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.5")
        else()
                add_definitions("--target=x86_64-apple-macos10.7")
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --target=x86_64-apple-macos10.7")
        endif()
        endif()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif (UNIX)
    message("Target OS     : Unix")
    # build as static binary to run on most machines
    if (NOT ALISTAT_FLAGS MATCHES "static")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
else()
    # Note that IQ-TREE has NOT been tested on other platforms
    message("Target OS     : Unknown and untested yet")
endif()


##################################################################
# Setup compiler, currently supported GCC, CLANG, MSVC, and ICC
##################################################################

set(GCC "FALSE")   #  GNU compiler
set(CLANG "FALSE") # Clang compiler
set(ICC "FALSE")   # Intel compiler
set(VCC "FALSE")   # MS Visual C Compiler, note that it is different from MSVC variable
set(CLANG_UNDER_VS "FALSE") #Clang compiler, used from inside Visual Studio
# using C++11 standard
# disable AVX for NEON
if (__ARM_NEON)
    set(ALISTAT_FLAGS "${ALISTAT_FLAGS} novx")
elseif (CMAKE_CXX_COMPILER MATCHES "VISUAL STUDIO")
    set(CLANG_UNDER_VS "TRUE")
    #it won't recognize the -std=c++17 parameter.
    #Todo: don't hard-code this; figure out some way it can be passed in (though ideally, not the whole shebang).
    include_directories("C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\Llvm\\lib\\clang\\10.0.0\\include")
elseif(ALISTAT_FLAGS MATCHES "cpp14")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
    message("Compiler      : GNU Compiler (gcc)")
    set(GCC "TRUE")
#    set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -fms-extensions -Wno-deprecated")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g -ffunction-sections -fdata-sections")
    set(CMAKE_C_FLAGS_RELEASE "-O2 -g -ffunction-sections -fdata-sections")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-dead_strip")
    else()
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,--gc-sections")
    endif()
    # require at least gcc ${GCC_MIN_VERSION}
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_MIN_VERSION)
        message(FATAL_ERROR "GCC version must be at least ${GCC_MIN_VERSION}!")
    endif()
    if (WIN32)
        # disable AVX on Windows due to memory alignment
        set(ALISTAT_FLAGS "${ALISTAT_FLAGS} novx")
        message("WARNING: AVX is disabled on Windows as GCC does not properly suport memory alignment")
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Compiler      : Clang")
    set(CLANG "TRUE")
#    set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -Wno-nested-anon-types")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffunction-sections -fdata-sections")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ffunction-sections -fdata-sections")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-dead_strip")
    else()
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,--gc-sections")
    endif()

    # use libc++ per default in MacOS
    if (APPLE)
        SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()

    #remove -rdynamic for Clang under Linux
    if (UNIX AND ALISTAT_FLAGS MATCHES "static")
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(VCC "TRUE")
    message("Compiler      : MS Visual C++ Compiler")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message("Compiler      : Intel C++ Compiler (icc)")
    set(ICC "TRUE")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
else()
    message("Compiler      : Unknown and untested yet")
endif()
message("Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

set(EXE_SUFFIX "")

if (MSVC)
    # MS Visual Studio environment
    message("Exporting MS Visual Studio projects...")
    if (CLANG_UNDER_VS)
        #see https://clang.llvm.org/docs/UsersManual.html#clang-cl
        #note .GX is how you say -fexceptions
        add_definitions(/D_UWIN)
        set(CMAKE_C_FLAGS_RELEASE "/O2 /GX")
        set(CMAKE_C_FLAGS_DEBUG   "/D_UWIN /GX")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /I${EIGEN3_INCLUDE_DIR}")
    else()
        add_definitions(/MP) # enable multi-processor compilation
    endif()
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_definitions(/Ot /Oi)
        if (VCC)
            add_definitions(/O2)
        elseif (ICC)
            add_definitions(/O3)
        endif()
    endif()
endif()

# enable link time optimization
if (ALISTAT_FLAGS MATCHES "lto")
    #if (CLANG)
    #    set(COMBINED_FLAGS "${COMBINED_FLAGS} -flto=thin")
    #else()
        set(COMBINED_FLAGS "${COMBINED_FLAGS} -flto")
    #endif()
endif()

##################################################################
# Setup compiler flags
##################################################################

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -fno-default-inline -fno-inline -O2 -fno-omit-frame-pointer -g")
set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -O2 -fno-omit-frame-pointer -g")

if(CLANG AND ALISTAT_FLAGS MATCHES "static")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -Wl,--allow-multiple-definition")
endif()

if (ALISTAT_FLAGS MATCHES "libcxx")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()


if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message("C flags       : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    message("CXX flags     : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("C flags       : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    message("CXX flags     : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Profile")
    message("C flags       : ${CMAKE_C_FLAGS_PROFILE} ")
    message("CXX flags     : ${CMAKE_CXX_FLAGS_PROFILE} ")
endif()

message("LINKER flags  : ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")


if (GCC)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
    set(CMAKE_CXX_FLAGS_MEM "-g -O1")
    set(CMAKE_C_FLAGS_MEM "-g -O1")
elseif (CLANG AND NOT CLANG_UNDER_VS)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
    set(CMAKE_CXX_FLAGS_MEM "-g -O1")
    set(CMAKE_C_FLAGS_MEM "-g -O1")
endif()

# Check if __ARM_NEON is defined
#include(cmake_utls/check_neon_exists.cmake)
set (__ARM_NEON "FALSE")
set (NEON 0)
if (APPLE)
    execute_process(
      COMMAND uname -m
      RESULT_VARIABLE result
      OUTPUT_VARIABLE OSX_NATIVE_ARCHITECTURE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(OSX_NATIVE_ARCHITECTURE STREQUAL "arm64")
        message("APPLE SILICON (M1/M2/M3) DETECTED!")
        set (__ARM_NEON "TRUE")
        set (NEON 1)
    endif()
elseif (UNIX AND NOT APPLE)  # Unix and Linux
    execute_process(
      COMMAND uname -m
      RESULT_VARIABLE result
      OUTPUT_VARIABLE OS_NATIVE_ARCHITECTURE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(OS_NATIVE_ARCHITECTURE STREQUAL "aarch64")
        message("Unix/Linux ARM64 DETECTED!")
        set (__ARM_NEON "TRUE")
        set (NEON 1)
    endif()
endif()

##################################################################
# setup linking flags
##################################################################

# link special lib for WIN32
if (WIN32)
    set(PLATFORM_LIB "ws2_32")
else()
    set(PLATFORM_LIB "m")
endif()

if (ALISTAT_FLAGS MATCHES "libcxx")
    set(STD_LIB "c++abi")
endif()

set(THREAD_LIB "")
if (NOT ALISTAT_FLAGS MATCHES "single")
    if (MSVC)
        if (BINARY32)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LIBPATH:${PROJECT_SOURCE_DIR}/lib32")
            set(THREAD_LIB "pthreadVC2")
        else()
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LIBPATH:${PROJECT_SOURCE_DIR}/lib")
            set(THREAD_LIB "pthreadVC2")
        endif()
    endif()

    if (CLANG AND BINARY32)
        set (ATOMIC_LIB "atomic")
    endif()

endif()


# basic linking libraries
target_link_libraries(alistat PUBLIC ${PLATFORM_LIB} ${STD_LIB} ${THREAD_LIB} ${ATOMIC_LIB})

# setup the executable name
##################################################################
set_target_properties(alistat PROPERTIES OUTPUT_NAME "alistat${EXE_SUFFIX}")

if (MSVC)
    set (BINARY_DIR "${PROJECT_BINARY_DIR}/Release")
else()
    set (BINARY_DIR "${PROJECT_BINARY_DIR}")
endif()

if (WIN32)
    if (MSVC)
        if (CLANG_UNDER_VS)
            ADD_CUSTOM_COMMAND(TARGET alistat POST_BUILD COMMAND copy "alistat${EXE_SUFFIX}.exe" "alistat${EXE_SUFFIX}-click.exe")
        else()
            ADD_CUSTOM_COMMAND(TARGET alistat POST_BUILD COMMAND copy "Release\\alistat${EXE_SUFFIX}.exe" "Release\\alistat${EXE_SUFFIX}-click.exe")
        endif()
    else()
        ADD_CUSTOM_COMMAND(TARGET alistat POST_BUILD COMMAND copy "alistat${EXE_SUFFIX}.exe" "alistat${EXE_SUFFIX}-click.exe")
    endif()
endif()

##############################################################
# add the install targets
##############################################################
install (TARGETS alistat DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/test.fas" DESTINATION examples)
install (FILES "${PROJECT_SOURCE_DIR}/M_Ento.fas" DESTINATION examples)
install (FILES "${PROJECT_SOURCE_DIR}/M_Ento.partition" DESTINATION examples)
install (FILES "${PROJECT_SOURCE_DIR}/AliStat_manual.pdf" DESTINATION .)

if (WIN32)
    install (FILES "${BINARY_DIR}/alistat${EXE_SUFFIX}-click.exe" DESTINATION bin)
    if (NOT ALISTAT_FLAGS MATCHES "single" AND MSVC)
        if (BINARY32)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib32/pthreadVC2.dll" DESTINATION bin)
        else()
            install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadVC2.dll" DESTINATION bin)
        endif()
    endif()
endif()

###############################################################################
## packaging ##################################################################
###############################################################################

include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

if(WIN32 OR APPLE)
  set(CPACK_GENERATOR "ZIP")
  set(CPACK_SOURCE_GENERATOR "ZIP")
else()
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR "TGZ")
endif()

set (SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    if (IQTREE_FLAGS MATCHES "oldmac")
        set (SYSTEM_NAME "macOS10.5")
    elseif (OSX_NATIVE_ARCHITECTURE STREQUAL "arm64")
        set (SYSTEM_NAME "macOS-arm")
    else()
        set (SYSTEM_NAME "macOS-intel")
    endif()
endif()

if (NOT APPLE AND UNIX)
    if (__ARM_NEON)
        set (SYSTEM_NAME "Linux-arm")
    else()
        set (SYSTEM_NAME "Linux-intel")
    endif()
endif()

if (BINARY32)
    set (SYSTEM_NAME "${SYSTEM_NAME}32")
endif()

if (ALISTAT_FLAGS MATCHES "KNL")
    set (SYSTEM_NAME "${SYSTEM_NAME}KNL")
endif()

set(CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}${EXE_SUFFIX}-${CPACK_PACKAGE_VERSION}-${SYSTEM_NAME}")

if (NOT APPLE)
    set(CPACK_STRIP_FILES TRUE)
endif()

include(CPack)
