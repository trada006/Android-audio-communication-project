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
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES=

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C \
	${TESTDIR}/TestFiles/f9 \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C \
	${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-std=c++0x -Ofast
CXXFLAGS=-std=c++0x -Ofast

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libc.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libc.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libc.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libc.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libc.a

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/BrentTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/DirectTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/FinitePulseTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/MasseyTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/NaiveTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/PilotOnlyTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/f9: ${TESTDIR}/tests/QAMTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f9 $^ ${LDLIBSOPTIONS} 

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/RecursiveTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C: ${TESTDIR}/tests/UtilTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc}   -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C $^ ${LDLIBSOPTIONS} 


${TESTDIR}/tests/BrentTest.o: tests/BrentTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/BrentTest.o tests/BrentTest.cpp


${TESTDIR}/tests/DirectTest.o: tests/DirectTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/DirectTest.o tests/DirectTest.cpp


${TESTDIR}/tests/FinitePulseTest.o: tests/FinitePulseTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/FinitePulseTest.o tests/FinitePulseTest.cpp


${TESTDIR}/tests/MasseyTest.o: tests/MasseyTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/MasseyTest.o tests/MasseyTest.cpp


${TESTDIR}/tests/NaiveTest.o: tests/NaiveTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/NaiveTest.o tests/NaiveTest.cpp


${TESTDIR}/tests/PilotOnlyTest.o: tests/PilotOnlyTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/PilotOnlyTest.o tests/PilotOnlyTest.cpp


${TESTDIR}/tests/QAMTest.o: tests/QAMTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -I. -I. -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/QAMTest.o tests/QAMTest.cpp


${TESTDIR}/tests/RecursiveTest.o: tests/RecursiveTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/RecursiveTest.o tests/RecursiveTest.cpp


${TESTDIR}/tests/UtilTest.o: tests/UtilTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/UtilTest.o tests/UtilTest.cpp


# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	    ${TESTDIR}/TestFiles/f9 || true; \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	    ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/C || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libc.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
