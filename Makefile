#
#  There exist several targets which are by default empty and which can be 
#  used for execution of your targets. These targets are usually executed 
#  before and after some main targets. They are: 
#
#     .build-pre:              called before 'build' target
#     .build-post:             called after 'build' target
#     .clean-pre:              called before 'clean' target
#     .clean-post:             called after 'clean' target
#     .clobber-pre:            called before 'clobber' target
#     .clobber-post:           called after 'clobber' target
#     .all-pre:                called before 'all' target
#     .all-post:               called after 'all' target
#     .help-pre:               called before 'help' target
#     .help-post:              called after 'help' target
#
#  Targets beginning with '.' are not intended to be called on their own.
#
#  Main targets can be executed directly, and they are:
#  
#     build                    build a specific configuration
#     clean                    remove built files from a configuration
#     clobber                  remove all built files
#     all                      build all configurations
#     help                     print help mesage
#  
#
#  Available make variables:
#
#     CND_BASEDIR                base directory for relative paths
#     CND_DISTDIR                default top distribution directory (build artifacts)
#     CND_BUILDDIR               default top build directory (object files, ...)
#     CONF                       name of current configuration
#     CND_PLATFORM_${CONF}       platform name (current configuration)
#     CND_ARTIFACT_DIR_${CONF}   directory of build artifact (current configuration)
#     CND_ARTIFACT_NAME_${CONF}  name of build artifact (current configuration)
#     CND_ARTIFACT_PATH_${CONF}  path to build artifact (current configuration)
#     CND_PACKAGE_DIR_${CONF}    directory of package (current configuration)
#     CND_PACKAGE_NAME_${CONF}   name of package (current configuration)
#     CND_PACKAGE_PATH_${CONF}   path to package (current configuration)
#
# NOCDDL

define newline


endef

# Environment 
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin

GFORT = gfortran -mcmodel=medium -fbacktrace
IFORT = ifort -r8 -mcmodel medium -traceback -diag-disable=10448
FORT = $(IFORT)
C = gcc
CPP = g++

BUILD_DIR = ./bin
LIB_DIR = ./lib
C_SRC_DIR = ./src/c
CPP_SRC_DIR = ./src/cpp
FORT_SRC_DIR = ./src/fortran
FORT_MICRO = micro_rates.f90
FORT_MACRO = macro_diffuse_into_and_along__internal \
             macro_diffuse_into_and_along__external \
            #  macro_Q2_diffuse_into \
            #  macro_Q2_diffuse_along \
            #  macro_Q2_always_rebind \
            #  macro_diffuse_into_and_along_slow_micro__external

#FILES =  ${FOLDER}macro_Q2.cpp  ${FOLDER}kiss.h

C_HEADERS = $(C_SRC_DIR)/all.h \
			$(C_SRC_DIR)/initializeData.h \
			$(C_SRC_DIR)/kiss.h \
			$(C_SRC_DIR)/nodeGrid.h \
			$(C_SRC_DIR)/parameters.h \
			$(C_SRC_DIR)/retrieve.h \
			$(C_SRC_DIR)/simulation.h \
			$(C_SRC_DIR)/test.h \
			$(C_SRC_DIR)/transfer.h \
			$(C_SRC_DIR)/timeList.h
C_SOURCE = $(C_SRC_DIR)/initializeData.c \
			$(C_SRC_DIR)/kiss.c \
			$(C_SRC_DIR)/main.c \
			$(C_SRC_DIR)/nodeGrid.c \
			$(C_SRC_DIR)/parameters.c \
			$(C_SRC_DIR)/retrieve.c \
			$(C_SRC_DIR)/simulation.c \
			$(C_SRC_DIR)/test.c \
			$(C_SRC_DIR)/transfer.c \
			$(C_SRC_DIR)/timeList.c
	
# build
build: .build-pre fort

fort: fort-micro fort-macro

c: c-macro

fort-micro: $(BUILD_DIR)/micro_rates # $(BUILD_DIR)/micro_wrapped

fort-macro: $(foreach file,$(FORT_MACRO),$(BUILD_DIR)/$(file))

f-macro-normal: $(BUILD_DIR)/macro-normal

f-macro-array: $(BUILD_DIR)/macro-array

shared: $(LIB_DIR)/kiss.so

$(BUILD_DIR)/micro_rates: $(FORT_SRC_DIR)/micro_rates.f90 $(BUILD_DIR)/kiss.o
	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/micro_rates.f90 -o $(BUILD_DIR)/micro_rates
    
# $(BUILD_DIR)/micro_wrapped: $(FORT_SRC_DIR)/micro_wrapped.f90 $(BUILD_DIR)/kiss.o
# 	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/micro_wrapped.f90 -o $(BUILD_DIR)/micro_wrapped

$(BUILD_DIR)/macro_diffuse_into_and_along__external: $(FORT_SRC_DIR)/macro_diffuse_into_and_along__external.f90 $(BUILD_DIR)/kiss.o
	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_diffuse_into_and_along__external.f90 -o $(BUILD_DIR)/macro_diffuse_into_and_along__external
    
# $(BUILD_DIR)/macro_diffuse_into_and_along_slow_micro__external: $(FORT_SRC_DIR)/macro_diffuse_into_and_along_slow_micro__external.f90 $(BUILD_DIR)/kiss.o
# 	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_diffuse_into_and_along_slow_micro__external.f90 -o $(BUILD_DIR)/macro_diffuse_into_and_along_slow_micro__external
    
$(BUILD_DIR)/macro_diffuse_into_and_along__internal: $(FORT_SRC_DIR)/macro_diffuse_into_and_along__internal.f90 $(BUILD_DIR)/kiss.o
	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_diffuse_into_and_along__internal.f90 -o $(BUILD_DIR)/macro_diffuse_into_and_along__internal
        
# $(BUILD_DIR)/macro_Q2_diffuse_into: $(FORT_SRC_DIR)/macro_Q2_diffuse_into.f90 $(BUILD_DIR)/kiss.o
# 	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_Q2_diffuse_into.f90 -o $(BUILD_DIR)/macro_Q2_diffuse_into
    
# $(BUILD_DIR)/macro_Q2_diffuse_along: $(FORT_SRC_DIR)/macro_Q2_diffuse_along.f90 $(BUILD_DIR)/kiss.o
# 	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_Q2_diffuse_along.f90 -o $(BUILD_DIR)/macro_Q2_diffuse_along
    
# $(BUILD_DIR)/macro_Q2_always_rebind: $(FORT_SRC_DIR)/macro_Q2_always_rebind.f90 $(BUILD_DIR)/kiss.o
# 	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_Q2_always_rebind.f90 -o $(BUILD_DIR)/macro_Q2_always_rebind
    
# $(BUILD_DIR)/macro-normal: $(FORT_SRC_DIR)/macro_brad_scratch.f90 $(BUILD_DIR)/kiss.o
# 	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_brad_scratch.f90 -o $(BUILD_DIR)/macro-normal
    
# $(BUILD_DIR)/macro-array: $(FORT_SRC_DIR)/macro_rng_array.f90 $(BUILD_DIR)/kiss.o
# 	$(FORT) $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/macro_rng_array.f90 -o $(BUILD_DIR)/macro-array
    
cpp: $(BUILD_DIR)/cpp_macro_Q2

$(BUILD_DIR)/cpp_macro_Q2: $(CPP_SRC_DIR)/macro_Q2.cpp $(BUILD_DIR)/kiss.o
	$(CPP) -std=c++11 -o $(BUILD_DIR)/cpp_macro_Q2 $(CPP_SRC_DIR)/macro_Q2.cpp $(BUILD_DIR)/kiss.o

$(LIB_DIR)/kiss.so: $(C_SRC_DIR)/kiss.c
	$(C) -fPIC -std=c99 -shared -o $(LIB_DIR)/kiss.so $(C_SRC_DIR)/kiss.c

$(BUILD_DIR)/kiss.o: $(C_SRC_DIR)/kiss.c
	$(C) -c $(C_SRC_DIR)/kiss.c -o $(BUILD_DIR)/kiss.o
#	gcc -c -std=c99 $(C_SRC_DIR)/kiss.c -o $(BUILD_DIR)/kiss.o

c-macro: $(C_HEADERS) $(C_SOURCE) # This line sets what files make looks at to determine if it needs to recompile.
	mpicc -O3 -lm -std=c11 -o $(BUILD_DIR)/c_macro ${C_SOURCE}

.build-pre:

.build-post:
# Add your post 'build' code here...


# clean
clean: .clean-post
	rm -r $(BUILD_DIR)
	mkdir $(BUILD_DIR)

.clean-pre:
# Add your pre 'clean' code here...

.clean-post:
# Add your post 'clean' code here...


# clobber
clobber: .clobber-post

.clobber-pre:
# Add your pre 'clobber' code here...

.clobber-post:
# Add your post 'clobber' code here...


# all
all: .all-post

.all-pre:
# Add your pre 'all' code here...

.all-post:
# Add your post 'all' code here...


# build tests
build-tests: .build-tests-post

.build-tests-pre:
# Add your pre 'build-tests' code here...

.build-tests-post:
# Add your post 'build-tests' code here...


# run tests
test: .test-post

.test-pre: build-tests
# Add your pre 'test' code here...

.test-post:
# Add your post 'test' code here...


# help
help: .help-post

.help-pre:
# Add your pre 'help' code here...

.help-post:
# Add your post 'help' code here...


