override HOME := ./..
include $(HOME)/Makefile 

% : %.cpp Makefile $(OBJS) 
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) $< -o $@ 