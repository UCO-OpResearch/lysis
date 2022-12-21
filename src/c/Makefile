#############################################################
# How to use:

# To use go into the directory with the files and type make.
# It will default to using this file as long as the name is Makefile.
# Make will output some text telling you what it has done.
# If make detects that a change has been made, it will recompile the program.
# If make doesn't detect any changes, it will say that it is up-to-date.
# website for GNU Make https://www.gnu.org/software/make/

# Note: This is simplified for easier use. Doesn't really matter since the compilation takes very little time.

#############################################################
# Files

# new TimeList
#HEADERS = all.h dynamicArray.h initializeData.h kiss.h nodeGrid.h parameters.h retrieve.h simulation.h test.h transfer.h
#SOURCE = dynamicArray.c initializeData.c kiss.c main.c nodeGrid.c parameters.c retrieve.c simulation.c test.c transfer.c

# has old TimeList
HEADERS = all.h initializeData.h kiss.h nodeGrid.h parameters.h retrieve.h simulation.h test.h transfer.h timeList.h
SOURCE = initializeData.c kiss.c main.c nodeGrid.c parameters.c retrieve.c simulation.c test.c transfer.c timeList.c

#############################################################
#Compiler and compiler options

# Switched to gcc 4.7.2 because the version of gdb debugger program on Buddy doesn't work with 4.9.2
# C_COMPILER = /act/openmpi-1.8/gcc-4.7.2/bin/mpicc
C_COMPILER = /home/Bryan/shared/mpich-3.4/bin/mpicc
OPENMPI_COMPILER = /home/Bryan/shared/openmpi-4.1.0/bin/mpicc

# Compilation Options
# -std=c11 sets the c language standard to 2011
# use -ggdb for debugging, take out -O optimization when using -ggdb
# use -O3 for optimization. Use for better performance when not debugging.
# use -g -rdynamic -D BACKTRACE when using the automatic backtrace.
	# -rdynamic loads additional symbols for backtrace
	# when an error occurs the list of functions is output to the slurm output file.
	# note that you can also do a backtrace in gdb
# -lm links the math.h library 
	# note that the -lm option needs to be before the source files. Otherwise the files that need the math library
	# will not be able to find the library.
# -D NEWTIMELIST for new TimeList

# These are the typical options for normal runs
OPTIONS = -O3 -lm -std=c11

# For typical debugging
DEBUG_OPTIONS = -ggdb -lm -std=c11

# For use with backtrace
BACKTRACE_OPTIONS = -g -rdynamic -lm -std=c11 -D BACKTRACE

#############################################################

ParallelMacro: ${HEADERS} ${SOURCE} # This line sets what files make looks at to determine if it needs to recompile.
	${C_COMPILER} ${OPTIONS} -o ParallelMacro ${SOURCE}

backtrace: ${HEADERS} ${SOURCE}
	${C_COMPILER} ${BACKTRACE_OPTIONS} -o ParallelMacro ${SOURCE}

debug: ${HEADERS} ${SOURCE}
	${C_COMPILER} ${DEBUG_OPTIONS} -o ParallelMacro ${SOURCE}

openmpi: ${HEADERS} ${SOURCE}
	${OPENMPI_COMPILER} ${OPTIONS} -o ParallelMacro ${SOURCE}