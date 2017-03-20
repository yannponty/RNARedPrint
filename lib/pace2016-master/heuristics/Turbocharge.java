package heuristics;

import graph.Graph;
import graph.TreeGraph;
import graph.TreeVertex;
import heuristics.greedy.Greedy;
import heuristics.greedy.GreedyCriteria;
import heuristics.greedy.GreedyDegree;
import heuristics.greedy.GreedyFillIn;

import java.util.*;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: heuristics
 *         Description: The turbocharge heuristic.
 */
public class Turbocharge {
  private long seed;
  private static final boolean DEBUG = false;
  private static final int L = 8;

  private TreeGraph bestTree;
  private TreeGraph minDegreeTree;
  private TreeGraph minFillTree;
  private TreeGraph singleNodeTree;

  public int maxBacktrack;

  public Turbocharge(int numVertices) {
    this(new Random().nextLong(), numVertices);
  }

  public Turbocharge(long seed, int numVertices) {
    this.seed = seed;
    this.maxBacktrack = -1;
    this.bestTree = null;
    this.minDegreeTree = null;
    this.minFillTree = null;

    // Construct a tree decomposition with a single vertex. This means every
    // vertex in the graph G is in a single bag. This tree is only ever used
    // if a SIGUSR1 or SIGTERM signal is sent before any reasonable tree
    // decomposition is constructed.
    this.singleNodeTree = new TreeGraph();
    TreeVertex singleNode = new TreeVertex("0", -1);
    for (int i = 0; i < numVertices; i++) {
      singleNode.addToBag(i);
    }
    this.singleNodeTree.addTreeVertex(singleNode);
  }

  public long getSeed() {
    return this.seed;
  }

  /**
   * @return Returns the best tree found so far during an execution
   * binaryTurboCharge.
   */
  public TreeGraph bestTreeSoFar() {
    if (this.bestTree != null) {
      return this.bestTree;
    } else if (this.minDegreeTree != null && this.minFillTree != null) {
      if (this.minDegreeTree.width() <= this.minFillTree.width()) {
        return this.minDegreeTree;
      }
      return this.minFillTree;
    } else if (this.minDegreeTree != null || this.minFillTree != null) {
      if (this.minDegreeTree != null) {
        return this.minDegreeTree;
      }
      return this.minFillTree;
    }
    return this.singleNodeTree;
  }

  public static void printStatusUpdate(TreeGraph T) {
    System.out.printf("c status %d %d\n",
        T.width() + 1,
        System.currentTimeMillis());
  }

  /**
   * Runs the turbocharge algorithm on a Graph, using a biased binary search
   * to find the best value of k. Makes a call to both greedy min degree and
   * greedy min fill 'NUM_GREEDY_CALLS` times, and at most `NUM_TURBO_CALLS`
   * calls to the turbo algorithm.
   * @param G
   * @return The best tree decomposition found.
   */
  public TreeGraph binaryTurbocharge(Graph G) {
    // The maximum number of times we call the turbo algorithm.
    int NUM_TURBO_CALLS = 4;
    int NUM_GREEDY_CALLS = 3;

    if (DEBUG)
      System.out.println("Seed = " + this.seed);

    Random random = new Random(this.seed);

    // Run the min degree and min fill algorithm.
    if (DEBUG)
      System.out.println("Running min degree " + NUM_GREEDY_CALLS + " time(s)");
    GreedyCriteria minDegree = new GreedyDegree(random, this.seed);
    int minDegreeTW = G.size() - 1;
    int[] minDegreePermutation = new int[G.size()];
    for (int i = 1; i <= NUM_GREEDY_CALLS; i++) {
      ArrayList<Integer> degreeOrdering =
          Greedy.greedyWithCriteria(G, minDegree);
      TreeGraph tree =
          ConstructTree.orderingToTreeGraph(G, degreeOrdering);
      int width = tree.width();
      if (width < minDegreeTW) {
        minDegreeTW = width;
        this.minDegreeTree = tree;
        System.arraycopy(G.permutation, 0,
            minDegreePermutation, 0, G.permutation.length);
        printStatusUpdate(this.minDegreeTree);
      }
      if (DEBUG) System.out.println("==> " + width);
    }


    if (DEBUG)
      System.out.println("Running min fill " + NUM_GREEDY_CALLS + " time(s)");
    GreedyCriteria minFill = new GreedyFillIn(random, this.seed);
    int minFillTW = G.size() - 1;
    int[] minFillPermutation = new int[G.size()];
    for (int i = 1; i <= NUM_GREEDY_CALLS; i++) {
      ArrayList<Integer> fillOrdering =
          Greedy.greedyWithCriteria(G, minFill);
      TreeGraph tree =
          ConstructTree.orderingToTreeGraph(G, fillOrdering);
      int width = tree.width();
      if (width < minFillTW) {
        minFillTW = width;
        this.minFillTree = tree;
        System.arraycopy(G.permutation, 0,
            minFillPermutation, 0, G.permutation.length);

        if (minFillTW < minDegreeTW) {
          printStatusUpdate(this.minFillTree);
        }
      }
      if (DEBUG) System.out.println("==> " + width);
    }

    // Determine the best algorithm to use, breaking ties in minDegrees
    // favour because it's faster.
    GreedyCriteria alg;
    int k;
    if (minDegreeTW <= minFillTW) {
      k = minDegreeTW;
      alg = minDegree;
      this.bestTree = this.minDegreeTree;
      System.arraycopy(minDegreePermutation, 0,
          G.permutation, 0, minDegreePermutation.length);
    } else {
      k = minFillTW;
      alg = minFill;
      this.bestTree = this.minFillTree;
      System.arraycopy(minFillPermutation, 0,
          G.permutation, 0, minFillPermutation.length);
    }

    double min = 0.94;
    double max = roundToHalf((k-1)/(double)k);
    int bestTW = k;

    // We cannot do better than this.
    if (bestTW == 1) {
      return this.bestTree;
    }

    // If the interval [min, max] is non-existent, try k - 1, k - 2, etc.
    // Otherwise, binary search on [min, max] increments of 0.005.
    if (max <= min) {
      for (int i = 1; i <= NUM_TURBO_CALLS; i++) {
        int newK = bestTW - 1;
        ArrayList<Integer> ordering =
            this.findOrdering(G, newK, this.L, alg, false);
        // If we found an ordering, great! Otherwise, we cannot hope to do
        // better than this.
        if (ordering != null) {
          TreeGraph T = ConstructTree.orderingToTreeGraph(G, ordering);
          if (T.width() < bestTW) {
            bestTW = T.width();
            this.bestTree = T;
            printStatusUpdate(this.bestTree);
          }
        } else {
          break;
        }
      }
    } else {
      // Keep track of the upper and lowerbound k values.
      int lbk;
      int ubk;
      int turboCalls = 0;
      while (min < max) {
        if (++turboCalls > NUM_TURBO_CALLS) {
          break;
        }

        // Make sure we don't call turbo again with the same value of k.
        double mid = roundToHalf((min + max) / 2);
        int midK = Math.max((int) Math.ceil(mid * k), 1);
        lbk = Math.max((int) Math.ceil(min * k), 1);
        ubk = Math.max((int) Math.ceil(max * k), 1);
        if (midK <= lbk) {
          midK++;
          if (midK >= ubk) {
            break;
          }
        }
        if (midK >= ubk) {
          midK--;
          if (midK <= lbk) {
            break;
          }
        }

        ArrayList<Integer> ordering =
            this.findOrdering(G, midK, this.L, alg, false);
        boolean decreaseTW = false;
        // If we found an ordering, check if it's better than what we have..
        if (ordering != null) {
          TreeGraph T = ConstructTree.orderingToTreeGraph(G, ordering);
          if (T.width() < bestTW) {
            bestTW = T.width();
            this.bestTree = T;
            decreaseTW = true;
            printStatusUpdate(this.bestTree);
          }
        }

        // If we found a tree, then let's try and find a better one by
        // decreasing k. Otherwise, raise k so that we can find something.
        if (decreaseTW) {
          max = mid;
        } else {
          min = mid;
        }
      }
    }

    if (DEBUG) System.out.println("Found tree width: " + bestTW);

    return this.bestTree;
  }

  public static double roundToHalf(double d) {
    return Math.round(d * 200) / 200.0;
  }

  /**
   * Finds an ordering using the turobcharge heuristic. Returns
   * ordering with treewidth at most k, otherwise null.
   * @param G The graph G.
   * @param k The upperbound on the treewidth.
   * @param l The parameter l (how far we should backtrack).
   * @param greedy The GreedyCriteria to use when selecting the next vertex.
   * @param generatePerm True if the turbo algorithm should generate a *new*
   *                     permutation of the vertices of G, using the random
   *                     object from greedy.
   * @return An ordering with treewidth at most k, otherwise null.
   */
  public ArrayList<Integer> findOrdering(Graph G, int k, int l,
                                         GreedyCriteria greedy,
                                         boolean generatePerm) {
    if (DEBUG)
      System.out.println("Running turbocharge with k = " + k + ", l = " + l);
    if (generatePerm) {
      G.generateRandomPermutation(greedy.getRnd());
    }
    Graph H = G.clone();
    ArrayList<Integer> ordering = new ArrayList<Integer>();
    for (int i = 1; i <= G.size(); i++) {
      int v = greedy.getNextVertex(H);
      // If adding the vertex does not make the treewidth greater than k, add
      // it to our ordering. Otherwise, backtrack l steps.
      if (H.degree(v) <= k) {
        ordering.add(v);
        ConstructTree.eliminateVertex(H, v);
      } else {
        if (DEBUG)
          System.out.println("==> Backtracking..");
        // If l > i - 1, then we can only backtrack up to i - 1 vertices.
        int backtrackLength = Math.min(l, i - 1);

        Graph backtrackG = G.clone();
        // Eliminate the first i - l - 1 choices.
        for (int j = 0; j < i - 1 - backtrackLength; j++) {
          ConstructTree.eliminateVertex(backtrackG, ordering.get(j));
        }
        ArrayList<Integer> W = new ArrayList<Integer>();
        for (int j = 0; j < backtrackG.size(); j++) {
          int u = backtrackG.permutation[j];
          if (backtrackG.contains(u) && backtrackG.degree(u) <= k) {
            W.add(u);
          }
        }
        ArrayList<Integer> extendedOrdering =
            ICTreeWidth(backtrackG, W, k, backtrackLength + 1,
                        greedy.getCriteriaAsString());
        if (extendedOrdering == null) {
          return null;
        }
        Collections.reverse(extendedOrdering);

        H = backtrackG.clone();
        // Keep a record of how much we actually backtrack.
        int backtracking = -1;
        boolean startOfNewOrdering = false;
        // Update the ordering to have a list of size
        // (i - l - 1) + (l + 1) = i, and the eliminate the new vertices
        // chosen in H.
        for (int j = i - 1 - backtrackLength; j < i; j++) {
          int vertex = extendedOrdering.get(j - i + 1 + backtrackLength);
          ConstructTree.eliminateVertex(H, vertex);
          if (j >= i - 1) {
            ordering.add(vertex);
          } else {
            if (ordering.get(j) != vertex && !startOfNewOrdering) {
              // If the two vertices in the original and new ordering are
              // different, then everything else after this index is a new
              // ordering (i.e. this is the maximum backtrack).
              startOfNewOrdering = true;
              backtracking = i - 1 - j;
            }
            ordering.set(j, vertex);
          }
        }

        if (DEBUG)
          System.out.println("====> Backtracked " + Math.max(backtracking, 0));
        if (backtracking > maxBacktrack)
          maxBacktrack = backtracking;
      }
    }

    return ordering;
  }

  /**
   * IC-Treewidth, used by the turbocharged min-degree algorithm (see notes).
   * @param G The graph.
   * @param W Vertex set.
   * @param k Parameter k.
   * @param l Parameter l.
   * @param criteria A criteria used to determine how vertices are sorted.
   *                 Must be either "degree" or "fillin".
   * @return A **reverse** partial ordering of length l on the vertices of W,
   *         with treewidth at most k. Returns null is the case of a failure.
   */
  private ArrayList<Integer> ICTreeWidth(Graph G, ArrayList<Integer> W,
                                        int k, int l, String criteria) {
    // Sort vertices in W by the given criteria.
    if (criteria.equalsIgnoreCase("fillin")) {
      Graph.computeFillIn(G, W);
    }
    G.sortByCriteria(W, criteria);

    for (int v : W) {
      if (l == 1) {
        ArrayList<Integer> order = new ArrayList<Integer>();
        order.add(v);
        return order;
      }

      // Eliminate v from G to obtain a new graph. We record which edges were
      // added during this process, so that we can easily add v back into G.
      // Also can be used to keep track of the vertices that now have
      // degree <= k in H.
      ArrayList<Integer> newW = new ArrayList<Integer>(W);
      Set<Integer[]> edgesToRemove = new HashSet<Integer[]>();
      ArrayList<Integer> neighbours = new ArrayList<Integer>(G.neighboursOf(v));
      for (int i : neighbours) {
        for (int j : neighbours) {
          if (i == j)
            continue;
          if (!G.isNeighbour(i, j)) {
            G.addEdge(i, j);
            Integer[] edge = new Integer[2];
            edge[0] = i; edge[1] = j;
            edgesToRemove.add(edge);
          }
        }
      }
      for (int u : neighbours) {
        G.removeNeighbour(u, v);
        if (G.degree(u) <= k && !newW.contains(u))
          newW.add(u);
        else if (G.degree(u) > k)
          newW.remove(Integer.valueOf(u));
      }
      G.removeVertex(v);
      newW.remove(Integer.valueOf(v));

      ArrayList<Integer> partialOrdering =
          ICTreeWidth(G, newW, k, l - 1, criteria);

      // Undo the changes made by G. Adding the vertex v back, remove the
      // edges that were added, and restore the neighbourhood of v.
      G.addVertex(v);
      for (int neighbour : G.neighboursOf(v)) {
        G.addEdge(v, neighbour);
      }
      for (Integer[] edge : edgesToRemove) {
        G.removeEdge(edge[0], edge[1]);
      }

      if (partialOrdering != null) {
        partialOrdering.add(v);
        return partialOrdering;
      }
    }

    return null;
  }
}
