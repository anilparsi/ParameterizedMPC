cmake_minimum_required(VERSION 2.6)

PROJECT(pdMPC_discrete)

file(GLOB PROJECT_SOURCE 
	"src/*.h"
	"src/*.cpp")

include_directories(src)
	
add_executable(pdMPC_discrete ${PROJECT_SOURCE})
#add_library(pdMPC_discrete_lib STATIC ${PROJECT_SOURCE} )