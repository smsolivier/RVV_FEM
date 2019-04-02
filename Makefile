# make sure there aren't trailing spaces for any directory names! 

# root directory. use override HOME := in makefiles in exec or test directories 
# using override allows for relative paths 
HOME = .

# CXX = g++
CXX = riscv64-unknown-elf-g++

# directories with source files in them 
UTILS = $(HOME)/utils
FEM = $(HOME)/fem
LINALG = $(HOME)/linalg 
GENERAL = $(HOME)/general 
MESH = $(HOME)/mesh 

# # where make searches for source files 
VPATH = $(UTILS) $(FEM) $(LINALG) $(GENERAL) $(MESH) 
vpath %.s $(UTILS) $(FEM) $(LINALG) $(GENERAL) $(MESH) 

VARS = -DTESTING -DUSE_WARNINGS
ifeq ($(CXX), g++) 
else 
	VARS += -DUSE_RISCV 
	ASM = $(wildcard $(HOME)/linalg/*.s)
	ASSEMBLER = $(ASM) -Wa,-march=rv64imafdcv
endif 

OPT = -O3 -ffast-math -DNDEBUG 
# OPT = -g 

# # where to look for header files 
CFLAGS = -std=c++17 -I$(UTILS) -I$(FEM) -I$(LINALG) -I$(GENERAL) -I$(MESH) -I$(HOME) 
CFLAGS += $(OPT) $(VARS) 

# store object files and dependency files 
OBJ = $(HOME)/obj
DEP = $(HOME)/dep

# get all file names for all .cpp files 
SRCFILES = $(notdir $(wildcard $(HOME)/mesh/*.cpp $(HOME)/fem/*.cpp \
	$(HOME)/linalg/*.cpp $(HOME)/general/*.cpp))

# convert to object files and dependency files 
OBJS = $(patsubst %.cpp,$(OBJ)/%.o,$(SRCFILES))
DEPS = $(patsubst $(OBJ)/%.o,$(DEP)/%.d, $(OBJS))
TESTS = $(basename $(wildcard $(HOME)/test/*.cpp))

$(OBJ)/%.o : %.cpp $(HOME)/Makefile
	mkdir -p $(OBJ); $(CXX) -c $(CFLAGS) $(LIBS) $< -o $@ 
	mkdir -p $(DEP)
	$(CXX) -MM $(CFLAGS) $(LIBS) $< | sed -e '1s@^@$(OBJ)\/@' > $*.d; mv $*.d $(DEP)

all : $(OBJS) 
tests : $(TESTS)

-include $(DEPS)

cleanall : 
	rm -rf $(OBJ); rm -rf $(DEP); rm -f $(TESTS) 

listsrc : 
	@echo $(SRCFILES) 
listobj : 
	@echo $(OBJS)
listdep :
	@echo $(DEPS) 
listflags : 
	@echo $(CFLAGS)
listtests:
	@echo $(TESTS) 
listasm:
	@echo $(ASM) 
.PHONY : docs
docs : 
	cd $(HOME)/docs; doxygen Doxyfile