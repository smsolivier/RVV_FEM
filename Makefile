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
BUILD = $(HOME)/build

# # where make searches for source files 
VPATH = $(UTILS) $(FEM) $(LINALG) $(GENERAL) $(MESH) 
vpath %.S $(UTILS) $(FEM) $(LINALG) $(GENERAL) $(MESH) 

VARS = -DTESTING
# VARS += -DUSE_WARNINGS 
VARS += -DNDEBUG
# VARS += -DCH_NTIMER 
ifeq ($(CXX), g++) 
	EXE = ./
	VARS += -DUSE_UNWIND
	LIBS += -lunwind -ldl 
else 
	VARS += -DUSE_RISCV 
	VARS += -march=rv64gc
	ASM = $(wildcard $(HOME)/linalg/*.S $(HOME)/fem/*.S)
	ASSEMBLER = $(ASM) -Wa,-march=rv64gcv
	EXE = spike --isa=rv64gcv --dc=128:4:64 pk 
endif 

OPT = -O3 -ffast-math
# OPT += -funroll-loops
# OPT += -g 

# # where to look for header files 
CFLAGS = -std=c++17 -I$(UTILS) -I$(FEM) -I$(LINALG) \
	-I$(GENERAL) -I$(MESH) -I$(HOME) 
CFLAGS += $(OPT) $(VARS)

# store object files and dependency files 
OBJ = $(BUILD)/obj
DEP = $(BUILD)/dep

# get all file names for all .cpp files 
SRCFILES = $(notdir $(wildcard $(HOME)/mesh/*.cpp $(HOME)/fem/*.cpp \
	$(HOME)/linalg/*.cpp $(HOME)/general/*.cpp $(HOME)/utils/*.cpp))

# convert to object files and dependency files 
OBJS = $(patsubst %.cpp,$(OBJ)/%.o,$(SRCFILES))
DEPS = $(patsubst $(OBJ)/%.o,$(DEP)/%.d, $(OBJS))
TESTS = $(basename $(wildcard $(HOME)/test/*.cpp))
TESTOUT = $(addsuffix .out,$(TESTS)) 
EXES = $(basename $(wildcard $(HOME)/exe/*.cpp)) 

$(OBJ)/%.o : %.cpp $(HOME)/Makefile
	mkdir -p $(OBJ)
	$(CXX) -c $(CFLAGS) $(LIBS) $< -o $@ 
	mkdir -p $(DEP)
	$(CXX) -MM $(CFLAGS) $(LIBS) $< | sed -e '1s@^@$(OBJ)\/@' > $*.d; mv $*.d $(DEP)

all : $(OBJS)
tests : $(TESTS) 
run : $(TESTOUT) 

%.S : %.cpp $(HOME)/Makefile
	$(CXX) -S $(CFLAGS) -fverbose-asm $< -o $@ 

%.out : %
	$(EXE)$<

-include $(DEPS)

clean: 
	rm -rf $(BUILD)
	rm -f $(TESTS) 
	rm -f $(TESTOUT) 
	rm -f $(EXES) 

listsrc : 
	@echo $(SRCFILES) 
listobj : 
	@echo $(OBJS)
listdep :
	@echo $(DEPS) 
listflags : 
	@echo $(CFLAGS)
listtests:
	@echo $(TESTOUT) 
listasm:
	@echo $(ASM) 
.PHONY : docs
docs : 
	cd $(HOME)/docs; doxygen Doxyfile
