# Set the project name
set(PROJ SCClient)

if (CMAKE_LINK2DTUWEC)
	set(PROJ "${PROJ}_DTUWEC")
endif()
if (CMAKE_LINK2ROSCO)
	set(PROJ "${PROJ}_ROSCO")
endif()
if (NOT CMAKE_LINK2MPI) # VC edit
	set(PROJ "${PROJ}_Standalone")
endif()


project(${PROJ} LANGUAGES Fortran)


message("Configuring Project: " ${PROJECT_NAME})
message("Root Directory is: " ${ROOT_DIR})

# set source code
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE SRC 
					  "${SRC_DIR}/YawAngleFromFile.f90"
					  "${SRC_DIR}/SCClientSubs.f90"
					  "${SRC_DIR}/SCClient_DISCON.f90"					 
	)
	
				   
message("Including the utils: " "${ROOT_DIR}/utils/CMakeLists.txt") 
include(${ROOT_DIR}/utils/CMakeLists.txt)

