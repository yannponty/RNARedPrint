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

all: RNARedPrint

RNARedPrint: $(OBJS) $(MAIN_SOURCE) $(TDLIB_OBJS)
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

clean: 
	rm -f $(PRODUCED) 

