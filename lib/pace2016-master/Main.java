import graph.*;
import heuristics.Turbocharge;

import java.lang.StringBuilder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class Main {
  public static void main(String[] args) {
    // java Main -s [seed] < [file path]
    final String SEED_FLAG = "-s";
    Long seed = null;
    if (args.length > 1 && args[0].equalsIgnoreCase(SEED_FLAG)) {
      seed = Long.parseLong(args[1]);
    }

    // Read the file from stdin.
    Graph G = null;
    Scanner sc = new Scanner(System.in);
    boolean foundStartOfFile = false;
    int n = 0;

    while (sc.hasNextLine()) {
      String line = sc.nextLine();
      String[] data = line.split("\\s+");
      if (line.startsWith("p") && !foundStartOfFile) {
        // Set the number of vertices needed.
        n = Integer.parseInt(data[2]);
        G = new Graph(n);
        for (int i = 0; i < n; i++) {
          G.addVertex(i, Integer.toString(i+1));
        }
        foundStartOfFile = true;
      } else if (!line.startsWith("c") && !data[0].equals(data[1])) {
        int u = Integer.parseInt(data[0]) - 1;
        int v = Integer.parseInt(data[1]) - 1;
        G.addEdge(u, v);
      }
    }

    final int numVertices = n;
    final Turbocharge t;
    if (seed != null)
      t = new Turbocharge(seed, numVertices);
    else
      t = new Turbocharge(numVertices);

    // Handle SIGTERM.
    Runtime.getRuntime().addShutdownHook(new Thread() {
      @Override
      public void run() {
        TreeGraph tree = t.bestTreeSoFar();
        printTreeGraph(numVertices, tree);
      }
    });

    // Run the turbocharge algorithm.
    TreeGraph tree = t.binaryTurbocharge(G);
  }

  public static void printTreeGraph(int numVertices, TreeGraph T) {
    // If there are less than n - 1 edges (where n is the size of the tree
    // decomposition), find all connected components in the forest and add
    // arbitrary edges between the components to form a fully connected tree.
    if (T.edges() < T.size() - 1) {
      ArrayList<ArrayList<TreeVertex>> components = T.findComponents();
      // Pick one component as the "root", add an arbitrary edge from each
      // other component to the root.
      TreeVertex root = components.get(0).get(0);
      for (int i = 1; i < components.size(); i++) {
        T.addEdge(components.get(i).get(0), root);
      }
    }

    HashMap<Integer, TreeVertex> indexToVertex =
        new HashMap<Integer, TreeVertex>();
    // Print out the size of the tree and graph.
    System.out.printf("s td %d %d %d\n",
        T.size(),
        T.width() + 1,
        numVertices);

    // Print out all of the bags.
    for (int i = 0; i < T.size(); i++) {
      TreeVertex t = T.getTreeVertices().get(i);
      // Use StringBuilder here, as may have long strings to print.
      StringBuilder line = new StringBuilder(String.valueOf(i + 1));
      for (int v : t.bag()) {
        line.append(" " + Integer.toString(v + 1));
      }
      indexToVertex.put(i + 1, t);
      System.out.printf("b %s\n", line.toString());
    }

    // Print out the edges.
    for (int i = 1; i <= T.size(); i++) {
      TreeVertex x = indexToVertex.get(i);
      for (int j = 1; j < i; j++) {
        TreeVertex y = indexToVertex.get(j);
        if (x.getNeighbours().contains(y)) {
          System.out.printf("%d %d\n", j, i);
        }
      }
    }
  }
}
