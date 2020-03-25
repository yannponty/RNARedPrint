MAIN_SOURCE = src/RNARedPrint.cpp  
SOURCES = src/DP.cpp  src/Nucleotide.cpp  src/RNAStructure.cpp  src/TreeDecomposition.cpp  src/Utils.cpp src/EnergyModels.cpp

OBJS = $(SOURCES:.cpp=.o)
EXEC = src/RNARedPrint

#COMPILER = $(CXX) -g -static-libgcc -static-libstdc++
COMPILER = $(CXX) -O2 -g -std=c++11 -DNDEBUG
JAVA_COMPILER = javac

TDLIB_SOURCES = lib/libtw/TD.java
TDLIB_OBJS = lib/libtw/TD.class
TDLIB_LIB = lib/libtw/libtw.jar

PRODUCED = $(OBJS) $(EXEC) $(TDLIB_OBJS)

SUBDIRS = "./lib/2016-pace-challenge-master/" "../flow-cutter-pace16-master" \
"../pace2016-master" "../pacechallenge-master"

# default installation prefix
PREFIX=/usr/local
ABS_PREFIX=$(realpath $(PREFIX))

all: $(EXEC)

$(EXEC): $(OBJS) $(MAIN_SOURCE) $(TDLIB_OBJS)
	$(COMPILER) $(OBJS) $(MAIN_SOURCE) -o $(EXEC)

%.o: %.cpp
	$(COMPILER) -c $< -o $@

%.class: %.java
	$(JAVA_COMPILER) -classpath $(TDLIB_LIB) $<

subsystem:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) all; \
	done

install: all
	install -d $(PREFIX)/bin
	install -d $(PREFIX)/share/RNARedPrint/lib
	cp -r lib $(PREFIX)/share/RNARedPrint
	install src/RNARedPrint $(PREFIX)/share/RNARedPrint/RNARedPrint
	install scripts/design-energyshift.py -D $(PREFIX)/bin
	install scripts/design-multistate.py -D $(PREFIX)/bin
	install --mode 644 scripts/RNARedPrintStructure.py -D $(PREFIX)/bin
	install --mode 644 scripts/RNARedPrintSampler.py -D $(PREFIX)/bin
	echo '#!/bin/sh\n$(ABS_PREFIX)/share/RNARedPrint/RNARedPrint $$* --prefix $(ABS_PREFIX)/share/RNARedPrint/' > $(PREFIX)/bin/RNARedPrint
	chmod 755 $(ABS_PREFIX)/bin/RNARedPrint

clean: 
	rm -f $(PRODUCED) 
