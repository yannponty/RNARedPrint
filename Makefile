MAIN_SOURCE = src/RNARedPrint.cpp  
SOURCES = src/DP.cpp  src/Nucleotide.cpp  src/RNAStructure.cpp  src/TreeDecomposition.cpp  src/Utils.cpp

OBJS = $(SOURCES:.cpp=.o)
EXEC = bin/RNARedPrint

COMPILER = g++ -g -static-libgcc -static-libstdc++
JAVA_COMPILER = javac

TDLIB_SOURCES = lib/libtw/TD.java
TDLIB_OBJS = lib/libtw/TD.class
TDLIB_LIB = lib/libtw/libtw.jar

PRODUCED = $(OBJS) $(EXEC) $(TDLIB_OBJS)

all: RNARedPrint

RNARedPrint: $(OBJS) $(MAIN_SOURCE) $(TDLIB_OBJS)
	$(COMPILER) $(OBJS) $(MAIN_SOURCE) -o $(EXEC)

%.o: %.cpp
	$(COMPILER) -c $< -o $@

%.class: %.java
	$(JAVA_COMPILER) -classpath $(TDLIB_LIB) $<

clean: 
	rm -f $(PRODUCED) 

