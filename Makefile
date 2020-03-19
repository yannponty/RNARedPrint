MAIN_SOURCE = src/RNARedPrint.cpp  
SOURCES = src/DP.cpp  src/Nucleotide.cpp  src/RNAStructure.cpp  src/TreeDecomposition.cpp  src/Utils.cpp src/EnergyModels.cpp

OBJS = $(SOURCES:.cpp=.o)
EXEC = bin/RNARedPrint

#COMPILER = g++ -g -static-libgcc -static-libstdc++
COMPILER = g++ -O2 -g -std=c++11 -DNDEBUG
JAVA_COMPILER = javac

TDLIB_SOURCES = lib/libtw/TD.java
TDLIB_OBJS = lib/libtw/TD.class
TDLIB_LIB = lib/libtw/libtw.jar

PRODUCED = $(OBJS) $(EXEC) $(TDLIB_OBJS)

SUBDIRS = "./lib/2016-pace-challenge-master/" "../flow-cutter-pace16-master" \
"../pace2016-master" "../pacechallenge-master"

# default installation prefix
PREFIX=/usr/local

all: $(EXEC)

$(EXEC): $(OBJS) $(MAIN_SOURCE) $(TDLIB_OBJS)
	mkdir -p bin
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
	install -d $(PREFIX)/share/RNARedPrint/lib
	cp -r lib $(PREFIX)/share/RNARedPrint
	install -d $(PREFIX)/bin/.bin
	install bin/RNARedPrint $(PREFIX)/bin/.bin/RNARedPrint
	install scripts/design-energyshift.py -D $(PREFIX)/bin
	install scripts/design-multistate.py -D $(PREFIX)/bin
	install --mode 644 scripts/RNARedPrintStructure.py -D $(PREFIX)/bin
	install --mode 644 scripts/RNARedPrintSampler.py -D $(PREFIX)/bin
	echo '#!/bin/sh\n$(PREFIX)/bin/.bin/RNARedPrint $$* --prefix $(PREFIX)/share/RNARedPrint/' > $(PREFIX)/bin/RNARedPrint
	chmod 755 $(PREFIX)/bin/RNARedPrint

clean: 
	rm -f $(PRODUCED) 
