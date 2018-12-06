cmake_minimum_required(VERSION 2.6)

PROJECT(pMPC)

file(GLOB PROJECT_SOURCE 
	"src/*.h"
	"src/*.cpp")

include_directories(src)
	
add_executable(pMPC ${PROJECT_SOURCE})
#add_library(pMPC STATIC ${PROJECT_SOURCE} )