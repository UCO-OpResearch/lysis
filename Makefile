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


# Environment 
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin

BUILD_DIR = ./bin
CPP_SRC_DIR = ./src/cpp
FORT_SRC_DIR = ./src/fortran
FORT_MICRO = micro_rates.f90
FORT_MACRO = macro_Q2_diffuse_into_and_along.f90
#FILES =  ${FOLDER}macro_Q2.cpp  ${FOLDER}kiss.h

	
# build
build: .build-pre cpp fort

fort: fort-micro $(BUILD_DIR)/macro

fort-micro: $(BUILD_DIR)/micro_rates

$(BUILD_DIR)/micro_rates: $(FORT_SRC_DIR)/micro_rates.f90 $(BUILD_DIR)/kiss.o
	ifort -r8 -mcmodel medium $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/micro_rates.f90 -o $(BUILD_DIR)/micro_rates

$(BUILD_DIR)/macro: $(FORT_SRC_DIR)/$(FORT_MACRO) $(BUILD_DIR)/kiss.o
	ifort -r8 -mcmodel medium $(BUILD_DIR)/kiss.o $(FORT_SRC_DIR)/$(FORT_MACRO) -o $(BUILD_DIR)/macro

cpp: $(BUILD_DIR)/cpp_macro_Q2

$(BUILD_DIR)/cpp_macro_Q2: $(CPP_SRC_DIR)/macro_Q2.cpp $(BUILD_DIR)/kiss.o
	g++ -std=c++11 -o ./bin/cpp_macro_Q2 $(CPP_SRC_DIR)/macro_Q2.cpp $(BUILD_DIR)/kiss.o

$(BUILD_DIR)/kiss.o: $(FORT_SRC_DIR)/kiss.c
	gcc -c $(FORT_SRC_DIR)/kiss.c -o $(BUILD_DIR)/kiss.o



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


