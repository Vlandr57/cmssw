# Short Makefile to compile "root macros"

CXX = g++
LINKER = g++ -o
Project = toyfit

scramtag = $(shell cd $$CMSSW_BASE; scram tool tag $(1) $(2))

ROOTLIBS := $(shell root-config --evelibs)
ROOTFLAGS := $(shell root-config --cflags)
ROOTINC := -I$(call scramtag, root_interface, INCLUDE)

ROOFITLIBS := -L$(call scramtag, roofitcore,LIBDIR ) -lRooFit -lRooFitCore -lMinuit
ROOFITFLAGS :=
ROOFITINC := -I$(call scramtag, roofitcore, INCLUDE)


USERLDFLAGS := $(ROOTLIBS) $(ROOFITLIBS)
USERCXXFLAGS := -g -Wall -O2
USERCXXFLAGS += $(ROOTFLAGS) $(ROOFITFLAGS)
USERINCFLAGS := $(ROOTINC) $(ROOFITINC)

main: $(Project).C
	$(CXX) $(USERINCFLAGS) $(USERCXXFLAGS) $(USERLDFLAGS) $(Project).C -o $(Project) $(USERLDFLAGS)




.PHONY: clean
	clean:   
	rm $(Project)

