CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated -fpermissive -std=c++11

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIBS 		= $(ROOTLIBS) $(ROOTGLIBS) $(WCSIMDIR)/libWCSimRoot.so
#LIBS 		= $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(WCSIMDIR)/libWCSimRoot.so

INC = $(WCSIMDIR)/include
SRC= $(WCSIMDIR)/src

CXXFLAGS += -I$(SRC) -I$(INC) 
#CXXFLAGS += -I$(LIBS) -I$(SRC) -I$(INC) 

TARGET= format_data_from_WCSim

all: $(TARGET)
format_data_from_WCSim: format_data_from_WCSim.o

%: %.o
	@echo "Now make $@"
	@echo $(LIBS)
# @$(CPP) -o $@ $< $(LIBS)
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS)
	@echo "..Compile done! "

%.o: %.C
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core.* *.log
