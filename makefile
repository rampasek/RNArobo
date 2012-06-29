CC      = g++
CFLAGS  = -O3 -Wall -Wno-sign-compare #-Wextra
LD      = $(CC)
LFLAGS  = -Wall

# Name of the final executable
EXE     = rnarobo
# Makes a list of the source (.cpp) files
SRCDIR = source/
SRCS := $(wildcard $(SRCDIR)*.cpp)
OBJS := ${SRCS:.cpp=.o}
DEPS := $(OBJS:.o=.d)


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
	@$(LD) $(LFLAGS) $(OBJS) -o $@
	
-include $(DEPS)

install: $(EXE)
	@echo Copying $(EXE) to $(BINDIR)/.
	cp $(EXE) $(BINDIR)/
   
clean:
	@echo Making clean.
	@rm -rf $(EXE) $(OBJS) $(DEPS) $(EXE_PROFILER_NAME) 

