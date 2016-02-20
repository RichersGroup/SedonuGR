# list of subdirectories
OBJDIR=.obj
SUBDIRS=$(shell find . -not -path '*/\.*' -type d)
$(shell mkdir -p $(SUBDIRS:%=$(OBJDIR)/%) >/dev/null)

# group together the external libraries that will be linked in
ALL_LIBS = $(GSL_LIBS) $(LUA_LDFLAGS) $(NULIB_LIBS) $(HDF5LIBS) $(LUA_LIBS)
ALL_INCS = $(GSL_INCL) $(LUA_INCL) $(OTHER_INCL) $(HDF5INCS) $(SUBDIRS:%=-I%)

# list all the executables
EXE_DIR=../exe
EXE_SRC_DIR=executables
EXE_SOURCES=$(shell find $(EXE_SRC_DIR)/*)
EXE_OBJECTS=$(EXE_SOURCES:%.cpp=$(OBJDIR)/%.o)
EXECUTABLES=$(EXE_SOURCES:$(EXE_SRC_DIR)/%.cpp=$(EXE_DIR)/%)

# all objects, excluding those with main definitions
OTHER_SOURCES = $(shell find . -type f -name '*.cpp' ! -path './$(EXE_SRC_DIR)/*')
OTHER_OBJECTS = $(OTHER_SOURCES:./%.cpp=$(OBJDIR)/%.o)

# list of dependencies
DEPS = $(OTHER_OBJECTS:.o=.d) $(EXE_OBJECTS:.o=.d)
DEPFLAGS = -MMD -MP -MF $(OBJDIR)/$*.Td

#============================================#

.PHONY: all clean

# link the executables
all: $(EXECUTABLES)
$(EXE_DIR)/%: $(OBJDIR)/$(EXE_SRC_DIR)/%.o $(OTHER_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ALL_LIBS)

# compile all the required objects and their dependencies
# creating temporary file so no corrupted dependencies on failed make
$(OBJDIR)/%.o: %.cpp
$(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< $(ALL_INCS) -IINCLUDES -o $@ 
	mv -f $(OBJDIR)/$*.Td $(OBJDIR)/$*.d

# include dependency files so recompile if header updates
$(OBJDIR)/%.d: ;
.PRECIOUS: $(OBJDIR)/%.d
-include $(DEPS)

# remove all compiled files
clean:
	rm -rf $(OBJDIR)
	rm -f $(TARGET)
	rm -f $(TESTS)

