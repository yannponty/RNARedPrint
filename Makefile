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

all: $(EXEC)

# $(TDLIB_OBJS)
$(EXEC): $(OBJS) $(MAIN_SOURCE)
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
	install -d $(PREFIX)/share/RNARedPrint/lib/treewidth-java
	cp -r lib/treewidth-java $(PREFIX)/share/RNARedPrint/lib
	install src/RNARedPrint $(PREFIX)/share/RNARedPrint/RNARedPrint
	install scripts/design-energyshift.py $(PREFIX)/bin
	install scripts/design-multistate.py $(PREFIX)/bin
	install scripts/RNARedPrint $(PREFIX)/bin
	install --mode 644 scripts/RNARedPrintStructure.py $(PREFIX)/bin
	install --mode 644 scripts/RNARedPrintSampler.py $(PREFIX)/bin

clean: 
	rm -f $(PRODUCED) 
