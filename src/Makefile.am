AUTOMAKE_OPTIONS = foreign

# these are flags sent to the C++ compiler
# -Wall : I don't know what this does
# -Werror: this will treat undefined and unused variables as errors
AM_CXXFLAGS = -Wall -Werror -std=c++0x

# This flag will build 'testexternals.C' as a test for the build environment
BUILT_SOURCES = \
  testexternals.C

# This is where we set the directories where we put our neccessary header files
AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libCrossCount.la

noinst_HEADERS = \
	dontIncludeMe.h

libCrossCount_la_LDFLAGS = \
  -L$(ROOTSYS)/lib `root-config --libs`

AM_LDFLAGS = \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib \
  -L$(libdir)

# These sources are compiled into the final library
libCrossCount_la_SOURCES = \
  CrossCount.C \
  SimpleFunction.C \
  DictOutput.cxx

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libCrossCount.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# This is where the dictionary file is generated
DictOutput.cxx: \
	CrossCount.h \
	SimpleFunction.h \
	CrossCountLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict.*
