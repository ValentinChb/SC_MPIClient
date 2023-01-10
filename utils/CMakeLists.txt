#This file is copied from the DTU utils repository


#Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

#Removing ZERO_CHECK project
set(CMAKE_SUPPRESS_REGENERATION true)

# Set platform variable
if (NOT CMAKE_SIZEOF_VOID_P)
    set (CMAKE_SIZEOF_VOID_P 4)
endif()
math(EXPR platform_bits "${CMAKE_SIZEOF_VOID_P} * 8")
set(PLATFORM "${CMAKE_SYSTEM_NAME}${platform_bits}")

if (NOT ROOT_DIR)
    # set root dir if not set in parent CMakeLists.txt
    set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()
# message("-----------CMAKE_Fortran_FLAGS_DEBUG (initvalue): " ${CMAKE_Fortran_FLAGS_DEBUG})
# message("-----------CMAKE_Fortran_FLAGS_RELEASE (initvalue): " ${CMAKE_Fortran_FLAGS_RELEASE})
if(NOT SRC_DIR)
  message("Source directories is missing, set the current directory as source directory!" )
  set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()


message("==========================================")
message("Project name: " ${PROJECT_NAME})
message("Platform: " ${PLATFORM})
message("Root dir: " ${ROOT_DIR})
message("Source directories:" ${SRC_DIR})
message("==========================================")

# VC edit: remove dependencies on low-level fortran libraries (other files in utils folder) for now
# if(CMAKE_Fortran_COMPILER)
#     # Adding common FORTRAN sources
    
#     # common version and output(logging) files
#     file(GLOB_RECURSE DEFAULTSRC  "${ROOT_DIR}/utils/version/BuildInfo.f90"
#                                   "${ROOT_DIR}/utils/version/version.rc"
#                                   "${ROOT_DIR}/utils/utils/string_utils.f90"
#                                   "${ROOT_DIR}/utils/io/logging.f90"
#                                   "${ROOT_DIR}/utils/io/path_utils.f90"
#                                   "${ROOT_DIR}/utils/io/dll_utils.f90"
#         )
#     list (APPEND SRC ${DEFAULTSRC})

#     if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
#         # GNU GCC specific files
#         file(GLOB_RECURSE COMPILERSRC   "${ROOT_DIR}/utils/portability/filesystem_tools_gfortran.f90"
#                                     "${ROOT_DIR}/utils/portability/support_compiler_gfrotran.f90"
#                                     "${ROOT_DIR}/utils/portability/dfwin_gfortran.f90" # interface for both linux and windows
#                            )  
#     elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
#         # Intel fortran specific files
#         file(GLOB_RECURSE COMPILERSRC   "${ROOT_DIR}/utils/portability/filesystem_tools_intel.f90"
#                                     "${ROOT_DIR}/utils/portability/support_compiler_intel.f90"
#                            )
#         if (${PLATFORM} STREQUAL "Linux64")
#             list(APPEND COMPILERSRC "${ROOT_DIR}/utils/portability/dfwin_linux.f90") # interface for both linux and windows (but breaks the compilation order in VS)
#         endif()
    
#     endif()
#     list (APPEND SRC ${COMPILERSRC})
# elseif(CMAKE_CXX_COMPILER)
#     # Adding common CPP sources
#     file(GLOB_RECURSE DEFAULTSRC    "${ROOT_DIR}/utils/version/BuildInfo.cpp"
#                                     "${ROOT_DIR}/utils/version/version.rc"
#                                     "${ROOT_DIR}/utils/portability/export.h"
#                                     )
#     list (APPEND SRC ${DEFAULTSRC})
# endif()

message("SOURCE FILE LISTS:  ")
list(LENGTH SRC SRCLEN)
message("Number of source files: " ${SRCLEN})
list (JOIN SRC "\n" txt)
# message (${txt}) # fails in some cases if number of source files too high
message("==========================================")
message("LIB MODULES:  " ${LIB_MOD})
message("==========================================")


if (NOT LIB)
    set(LIB dll)
endif()

# Creating exe, or shared or static library target
if (NOT BUILD_TYPE)
    message (FATAL_ERROR "BUILD_TYPE must be specified")
elseif (${BUILD_TYPE} STREQUAL "EXE")
    add_executable (${PROJECT_NAME} ${SRC})
    target_compile_definitions(${PROJECT_NAME} PRIVATE 
        MAIN
    )
elseif((${BUILD_TYPE} STREQUAL "SHARED") OR (${BUILD_TYPE} STREQUAL "STATIC"))
    #Creating static library target
    add_library (${PROJECT_NAME} ${BUILD_TYPE} ${SRC})
else()
    message (FATAL_ERROR "BUILD_TYPE must be one of EXE, SHARED, STATIC")
endif()


# set the current source directory as an include directory
include_directories(${SRC_DIR})
if(CLUSTER)
	set(CLUSTER_ID	"_cluster")
endif()
if (${PLATFORM} STREQUAL "Linux64")
    set(os_dir linux)
    set(CMAKE_SHARED_LIBRARY_PREFIX "")
    set(CMAKE_EXECUTABLE_SUFFIX ".exe")
	set(CMAKE_SHARED_LIBRARY_SUFFIX "${CLUSTER_ID}.so")
    add_definitions(-D_LINUX)
    set(WARN " -warn ")
    set(CHECK " -check ")
else()
    set(WARN " /warn:")
    set(CHECK " /check:")
    if (${PLATFORM} STREQUAL "Windows64")
        set (CMAKE_SHARED_LIBRARY_PREFIX "")
		set (CMAKE_SHARED_LIBRARY_SUFFIX "${CLUSTER_ID}_64.dll")
        set(os_dir x64)
        add_definitions(-D_WIN64)
	else() # Windows32
        set (CMAKE_SHARED_LIBRARY_PREFIX "")
		set(CMAKE_SHARED_LIBRARY_SUFFIX "${CLUSTER_ID}.dll")
        add_definitions(-D_WIN32)
        set(os_dir win32)       
    endif()
endif()


if(CMAKE_Fortran_COMPILER)
    
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        # message ("Using GCC" )
        target_compile_options(${PROJECT_NAME} PRIVATE
            -cpp;
            -ffree-form;
            -ffree-line-length-none;
            -fPIC;
            -fcray-pointer;
        )
		if (${PLATFORM} STREQUAL "Linux64")
			target_compile_definitions(${PROJECT_NAME} PRIVATE 
			  _LINUX
			)
		else()
			target_compile_definitions(${PROJECT_NAME} PRIVATE 
			  _WIN32
			)
		endif()

        set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -fcheck=bounds")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -static -fcheck=bounds")
		
        #VC edit: add flag to fix error "invalid register for .seh_savexmm"
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-asynchronous-unwind-tables")
    
	elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        # message ("Using Intel fortran")

        target_compile_options(${PROJECT_NAME} PRIVATE
            # if windows then <windows arg, e.g. ",/MP"> else <linux arg, e.g. "-multiple-processes">
			# Creates multiple processes that can be used to compile large numbers of source files at the same time.
			# https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-multiple-processes-mp
			$<IF:$<PLATFORM_ID:Windows>,/MP,-multiple-processes>;
			# Specifies a multi-threaded, static run-time library.
			# https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-libs
			$<IF:$<PLATFORM_ID:Windows>,/libs:${LIB};,>
			# Tells the linker to search for unresolved references in a multithreaded run-time library.
			# https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-threads
			$<IF:$<PLATFORM_ID:Windows>,/,->threads;
            $<IF:$<PLATFORM_ID:Windows>,/,->fpp;
            $<IF:$<PLATFORM_ID:Windows>,/Q,->fp-stack-check;
            $<IF:$<PLATFORM_ID:Windows>,/fp:source,-fp-model source>;      # floating point model
            $<IF:$<PLATFORM_ID:Windows>,/,->traceback;    
			# Ignore Remark # 7841: DLL IMPORT/EXPORT is not supported on this platform
			$<IF:$<PLATFORM_ID:Windows>,/Qdiag-disable:7841,-diag-disable=7841>;
			# Ignore Warning # 4217:  locally defined symbol <name> imported in function <fname>
			$<IF:$<PLATFORM_ID:Windows>,/Qdiag-disable:4217,-diag-disable=4217>;
			# Ignore Warning # 5462:  Global name too long
			$<IF:$<PLATFORM_ID:Windows>,,-diag-disable=5462>;
            $<IF:$<PLATFORM_ID:Windows>,${CODECOV},>
            $<IF:$<PLATFORM_ID:Windows>,/Qmkl:sequential;,>           #Tells the compiler to link using the parallel libraries in the Intel® MKL.  https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-mkl-qmkl    
            $<IF:$<PLATFORM_ID:Windows>,,-fPIC;>
        )

        
        
        

        target_link_options(${PROJECT_NAME} PRIVATE
            #linker
			# Ignore Warning # 4217:  locally defined symbol <name> imported in function <fname>
      		#$<IF:$<PLATFORM_ID:Windows>,,-ignore=4217>
      		# LINK : warning LNK4098: defaultlib 'MSVCRT' conflicts with use of other libs; use /NODEFAULTLIB:library
      		$<IF:$<PLATFORM_ID:Windows>,/ignore:4098,>;
      		# Tells the compiler to link using the sequential libraries in the Intel® MKL.
      		# https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-mkl-qmkl            
			#Tells the compiler to link using the sequential libraries in the Intel® MKL.
			#https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-mkl-qmkl			
			$<IF:$<PLATFORM_ID:Windows>,,-mkl:sequential;>           
            $<IF:$<PLATFORM_ID:Windows>,/INCREMENTAL:NO;,>
        )


        if (NOT defaults_set)
            set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${WARN}all${CHECK}all ${CHECK}noarg_temp_created")
            set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}${WARN}interfaces${CHECK}pointer${CHECK}bound${CHECK}uninit")
        endif()
    endif()
    if (NOT defaults_set)
            
        if (${PLATFORM} STREQUAL "Linux64")
          SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,'$ORIGIN'/lib -Wl,-rpath -Wl,'$ORIGIN'/")
          SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath -Wl,'$ORIGIN'/lib -Wl,-rpath -Wl,'$ORIGIN'/")
        else()
          # Default Incremental disables traceback, so delete default /Incremental it
            string(REPLACE "/INCREMENTAL" "" CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
            string(REPLACE "/INCREMENTAL" "" CMAKE_MODULE_LINKER_FLAGS_DEBUG "${CMAKE_MODULE_LINKER_FLAGS_DEBUG}")
            string(REPLACE "/INCREMENTAL" "" CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")
        endif()     
    endif()


	
	# VC edit: add define variables for MPI, DTUWEC and ROSCO
	if (CMAKE_LINK2MPI)
		 SET(CMAKE_Fortran_FLAGS_DEBUG "-DMPI ${CMAKE_Fortran_FLAGS_DEBUG}")
		 SET(CMAKE_Fortran_FLAGS_RELEASE "-DMPI ${CMAKE_Fortran_FLAGS_RELEASE}")
	endif()
    if (CMAKE_LINK2DTUWEC)
        SET(CMAKE_Fortran_FLAGS_DEBUG "-DDTUWEC ${CMAKE_Fortran_FLAGS_DEBUG}")
        SET(CMAKE_Fortran_FLAGS_RELEASE "-DDTUWEC ${CMAKE_Fortran_FLAGS_RELEASE}")
    endif()	
	if (CMAKE_LINK2ROSCO)
		 SET(CMAKE_Fortran_FLAGS_DEBUG "-DROSCO ${CMAKE_Fortran_FLAGS_DEBUG}")
		 SET(CMAKE_Fortran_FLAGS_RELEASE "-DROSCO ${CMAKE_Fortran_FLAGS_RELEASE}")
	endif()	
    

    message("Fortran: " ${CMAKE_Fortran_COMPILER_ID})
    message("Fortran compiler: " ${CMAKE_Fortran_COMPILER})
    message("FLAGS_DEGUG: " ${CMAKE_Fortran_FLAGS_DEBUG})
    message("FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE})
    message("FLAGS_EXE_LINKER: " ${CMAKE_EXE_LINKER_FLAGS})
    message("CMAKE_EXE_LINKER_FLAGS_DEBUG: " ${CMAKE_EXE_LINKER_FLAGS_DEBUG})
    message("CMAKE_EXE_LINKER_FLAGS_RELEASE: " ${CMAKE_EXE_LINKER_FLAGS_RELEASE})
    message("==========================================")

elseif(CMAKE_CXX_COMPILER)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
       SET(ENABLE_CXX11 "-std=c++11")

       EXECUTE_PROCESS(COMMAND "${CMAKE_CXX_COMPILER} -dumpversion" OUTPUT_VARIABLE GCC_VERSION)
       if (GCC_VERSION VERSION_LESS 4.7)
          SET(ENABLE_CXX11 "-std=c++0x")
       endif()

       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ENABLE_CXX11}")
    else()
        add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    endif()
    
    message("CXX " ${CMAKE_CXX_COMPILER_ID})
    message("CXX compiler: " ${CMAKE_CXX_COMPILER})
    message("FLAGS: " ${CMAKE_CXX_FLAGS})
    message("==========================================")
endif()


if(CLUSTER)
  target_compile_definitions(${PROJECT_NAME} PRIVATE CLUSTER)
endif()

# Specify include directories to use when compiling a given target
target_include_directories(${PROJECT_NAME} PRIVATE
          ${PROJECT_BINARY_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${ROOT_DIR}/utils/portability
)

#VC edit: link to mpi library 
if (CMAKE_LINK2MPI)
	target_link_libraries(${PROJECT_NAME} ${ROOT_DIR}/custom-build/src/libmsmpi.a)
	message("Linking to MPI library")	
endif()

#VC edit: link to dtuwec library 
if (CMAKE_LINK2DTUWEC)
	target_link_libraries(${PROJECT_NAME} ${ROOT_DIR}/custom-build/src/libDTUWEC4SC.a)
	message("Linking to DTUWEC library")				
endif()

#VC edit: link to rosco library 
if (CMAKE_LINK2ROSCO)
	target_link_libraries(${PROJECT_NAME} ${ROOT_DIR}/custom-build/src/libROSCO.a)
	message("Linking to ROSCO library")				
endif()

# VC edit: switch off version control for now
# if ((${CMAKE_GENERATOR} MATCHES "Visual Studio") AND (NOT TARGET NUM_LIB))
#     add_custom_command(TARGET ${PROJECT_NAME}
#        PRE_BUILD
#        COMMAND ${ROOT_DIR}\\utils\\version\\version.bat $(platform) $(configuration)
#     )
# else()
#     message("Project Build folder: " ${PROJECT_BINARY_DIR})
#     if (${PLATFORM} STREQUAL "Linux64")
#         execute_process(COMMAND bash ${ROOT_DIR}/utils/version/version.sh ${PLATFORM} ${CMAKE_BUILD_TYPE}
#                         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
#                         )
#     else()
#         execute_process(COMMAND ${ROOT_DIR}/utils/version/version.bat ${PLATFORM} ${CMAKE_BUILD_TYPE}
#                         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
#                         )
#     endif()
# endif()

function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
        message("There is no target named '${tgt}'")
        return()
    endif()

    # this list of properties can be extended as needed
    set(CMAKE_PROPERTY_LIST SOURCE_DIR BINARY_DIR COMPILE_DEFINITIONS
             COMPILE_OPTIONS INCLUDE_DIRECTORIES LINK_LIBRARIES)

    message("Configuration for target ${tgt}")

    foreach (prop ${CMAKE_PROPERTY_LIST})
        get_property(propval TARGET ${tgt} PROPERTY ${prop} SET)
        if (propval)
            get_target_property(propval ${tgt} ${prop})
            message (STATUS "${prop} = ${propval}")
        endif()
    endforeach(prop)

endfunction(print_target_properties)


print_target_properties(${PROJECT_NAME})
set(defaults_set 1)
message (${CMAKE_MAKE_PROGRAM})
message("==========================================")
