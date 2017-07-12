#


# Environment
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin
#VPATH=Build/

CPP_DIR =  C++\ code/
CPP_SRC = macro_Q2_cpp.o kiss.o -lncurses
CPP_FLAGS = -std=c++11
#FILES =  ${FOLDER}macro_Q2.cpp  ${FOLDER}kiss.h

C_FLAGS = -w

FORTRAN_DIR = Fortran\ code/
NCURSES_DIR = ${FORTRAN_DIR}pdsrc/
FORTRAN_SRC = macro_Q2_fort.o kiss.o macros.o ncurses.o -lncurses
FORTRAN_FLAGS = -I${NCURSES_DIR} -J${NCURSES_DIR} \
		-ffree-line-length-none


# build
build: 		.pre-build build_c build_fortran

full:		clean build_c build_fortran .post-build

#${C_SRC}: %.o : %.c
#	gcc -c $< -o $@


macro_Q2_cpp.o:
	g++ ${CPP_FLAGS} -c ${CPP_DIR}macro_Q2.cpp -o ${VPATH}$@

kiss.o:
	gcc ${C_FLAGS} -c ${CPP_DIR}kiss.c -o ${VPATH}$@

build_c:	${CPP_SRC}
	g++ ${CPP_FLAGS} ${CPP_SRC} -o macro_Q2_cpp

macros.o:
	gcc ${C_FLAGS} -c ${NCURSES_DIR}macros.c -o ${VPATH}$@

ncurses.o:
	gfortran ${FORTRAN_FLAGS} -c ${NCURSES_DIR}ncurses.f90 -o ${VPATH}$@

macro_Q2_fort.o:
	gfortran ${FORTRAN_FLAGS} -c ${FORTRAN_DIR}macro_Q2.f90 -o ${VPATH}$@

build_fortran:	${FORTRAN_SRC}
	gfortran ${FORTRAN_FLAGS} ${FORTRAN_SRC} -o macro_Q2_fort


clean:
	rm -f *.o
	rm -f macro_Q2*

.pre-build:
	@if [ -e macro_Q2_fort.o ]; then \
		if [ "$(shell stat -c %Y macro_Q2_fort.o)" -lt "$(shell stat -c %Y ${FORTRAN_DIR}macro_Q2.f90)" ]; then \
			rm macro_Q2_fort.o; echo "Removing out-of-date Fortran source."; \
		fi; \
	fi
	@if [ "$(shell stat -c %Y macro_Q2_cpp.o)" -lt "$(shell stat -c %Y ${CPP_DIR}macro_Q2.cpp)" ]; then \
		rm macro_Q2_cpp.o; echo "Removing out-of-date C++ source."; \
	fi

.post-build:

