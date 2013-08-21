CC      = g++
CFLAGS  = -O3 -Wall -Wno-sign-compare -msse -msse2 #-Wextra
LD      = $(CC)
LFLAGS  = -Wall
SLIBS   =

# To get POSIX timers we need to link librt
UNAME_S = $(shell uname -s)
ifneq ($(UNAME_S),Darwin)
    SLIBS += -lrt
endif

# Name of the final executable
EXE     = rnarobo
# Makes a list of the source (.cpp) files
SRCDIR = source/
SRCS := $(wildcard $(SRCDIR)*.cpp)
OBJS := ${SRCS:.cpp=.o}
DEPS := $(OBJS:.o=.d)

GOOGLE_MALLOCLIB =  #libtcmalloc_minimal.a
EXE_PROFILER_NAME = rnarobo_prof

# Where the binary will be installed
BINDIR = /usr/local/bin


# default target
all: $(EXE)

# define a suffix rule for .cpp -> .o
%.o: %.cpp
	@echo Creating object file for $*...
	@$(CC) $(CFLAGS) -MD -o $@ -c $<

$(EXE): $(OBJS)
	@echo Linking $@.
	@$(LD) $(LFLAGS) $(OBJS) $(SLIBS) -o $@
	
-include $(DEPS)

install: $(EXE)
	@echo Copying $(EXE) to $(BINDIR)/.
	cp $(EXE) $(BINDIR)/
   
clean:
	@echo Making clean.
	@rm -rf $(EXE) $(OBJS) $(DEPS) $(EXE_PROFILER_NAME) 

               
profiler: $(OBJS) $(GOOGLE_MALLOCLIB)
	$(LD) $(OBJS) $(GOOGLE_MALLOCLIB) -o $(EXE_PROFILER_NAME) -lprofiler
	

# "testing" machinery
TESTDIR = tests/
CMP = $(TESTDIR)tester.pl
TEST_LOG = $(TESTDIR)test_`date +"%d%b%T"`.txt
SUBDIR = 
INPUT_DIR = $(TESTDIR)in/$(SUBDIR)
OUTPUT_DIR = $(TESTDIR)out/$(SUBDIR)
KEY_DIR = $(TESTDIR)key/$(SUBDIR)
DESC = $(wildcard $(INPUT_DIR)*.des)
SEQS = $(wildcard $(INPUT_DIR)*.fa*)
TFLAGS = -c

test: $(EXE)
	@echo Running $(EXE) with input from $(INPUT_DIR).
	@num='0';
	@log_file="$(TEST_LOG)"; \
	for des in $(DESC); do \
		for seq in $(SEQS);	do \
		  num=`expr $$num + 1`; \
		  input="$$des $$seq";  \
	    output="$(OUTPUT_DIR)output$$num.out";  \
	    key="$(KEY_DIR)output$$num.out";  \
			echo "test $$num:\n  IN: $$input\n  OUT: $$output"; \
	    /usr/bin/time -p ./$(EXE) $(TFLAGS) $$input > $$output 2>> $$output; \
	    if [ $$? -eq 0 ]; then \
         echo '  FINISHED'; \
         ./$(CMP) $$key $$output -a $$log_file ; \
         if [ $$? -eq 1 ]; then \
            echo '  TEST OK'; \
         else \
            echo '  TEST FAILED'; \
         fi \
      else \
         echo '  FAILED!'; \
   	  fi \
		done \
	done
	
