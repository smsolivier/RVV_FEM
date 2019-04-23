override HOME := ./..
include $(HOME)/Makefile 

% : %.cpp Makefile $(OBJS) $(ASM) 
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) $(ASSEMBLER) $< -o $@ 
