CC= g++
INCLDIR= ./include/
SRCDIR= ./src/
OBJDIR= ./objs/
SRIMDIR = ./srim/
ROOT= `root-config --cflags --glibs`
CFLAGS= -g -Wall $(ROOT)
CPPFLAGS= -I$(INCLDIR) -I./ 
LDFLAGS= -L$(INCLDIR) $(ROOT)
SRC= $(wildcard $(SRCDIR)*.cpp)
OBJS= $(SRC:$(SRCDIR)%.cpp=$(OBJDIR)%.o)
DICT_PAGES=$(INCLDIR)nucleus.h $(INCLDIR)LinkDef_sim.h
DICT=$(SRCDIR)sim_dict.cxx
LIB= $(OBJDIR)sim_dict.o
EXE= simulator
SRIM_EXE= srim_clean


.PHONY: all clean

all: $(EXE) $(SRIM_EXE)

$(EXE): $(LIB) $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@

$(LIB): $(DICT)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $^
	mv $(SRCDIR)*.pcm ./

$(DICT): $(DICT_PAGES)
	rootcint -f $@ -c $^

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ 

$(SRIM_EXE): $(SRIMDIR)SRIMscript.cpp
	$(CC) -o $@ $^

clean:
	$(RM) $(OBJS) $(LIB) $(EXE) $(DICT) *.pcm
