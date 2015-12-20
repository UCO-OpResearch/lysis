#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/C++\ code/kiss.o \
	${OBJECTDIR}/C++\ code/macro_Q2.o \
	${OBJECTDIR}/C++\ code/tutorial.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/opt/local/include/boost -lboost_atomic-mt -lboost_chrono-mt -lboost_chrono-mt -lboost_container-mt -lboost_context-mt -lboost_coroutine-mt -lboost_date_time-mt -lboost_exception-mt -lboost_filesystem-mt -lboost_graph-mt -lboost_iostreams-mt -lboost_locale-mt -lboost_log-mt -lboost_log_setup-mt -lboost_math_c99-mt -lboost_math_c99f-mt -lboost_math_c99l-mt -lboost_math_tr1-mt -lboost_math_tr1f-mt -lboost_math_tr1l-mt -lboost_prg_exec_monitor-mt -lboost_program_options-mt -lboost_python-mt -lboost_random-mt -lboost_regex-mt -lboost_serialization-mt -lboost_signals-mt -lboost_system-mt -lboost_system-mt -lboost_test_exec_monitor-mt -lboost_thread-mt -lboost_timer-mt -lboost_timer-mt -lboost_unit_test_framework-mt -lboost_wave-mt -lboost_wserialization-mt -l%libname%

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/bloodclotting

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/bloodclotting: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/bloodclotting ${OBJECTFILES} ${LDLIBSOPTIONS}

.NO_PARALLEL:${OBJECTDIR}/C++\ code/kiss.o
${OBJECTDIR}/C++\ code/kiss.o: C++\ code/kiss.c 
	${MKDIR} -p ${OBJECTDIR} code
	${RM} "$@.d"
	$(COMPILE.c) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/C++\ code/kiss.o C++\ code/kiss.c

.NO_PARALLEL:${OBJECTDIR}/C++\ code/macro_Q2.o
${OBJECTDIR}/C++\ code/macro_Q2.o: C++\ code/macro_Q2.cpp 
	${MKDIR} -p ${OBJECTDIR} code
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/local/include/boost -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/C++\ code/macro_Q2.o C++\ code/macro_Q2.cpp

.NO_PARALLEL:${OBJECTDIR}/C++\ code/tutorial.o
${OBJECTDIR}/C++\ code/tutorial.o: C++\ code/tutorial.cpp 
	${MKDIR} -p ${OBJECTDIR} code
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/local/include/boost -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/C++\ code/tutorial.o C++\ code/tutorial.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/bloodclotting

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
