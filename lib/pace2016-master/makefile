JFLAGS =
JC = javac
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) $*.java 2> /dev/null

CLASSES = \
	graph/Graph.java \
	graph/TreeVertex.java \
	graph/TreeGraph.java \
	heuristics/ConstructTree.java \
	heuristics/greedy/GreedyCriteria.java \
	heuristics/greedy/GreedyDegree.java \
	heuristics/greedy/GreedyFillIn.java \
	heuristics/greedy/Greedy.java \
	heuristics/Turbocharge.java \
	Main.java \

all: classes

classes: $(CLASSES:.java=.class)

clean:
	$(RM) *.class graph/*.class \
	heuristics/*.class heuristics/greedy/*.class
