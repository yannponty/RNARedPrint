package heuristics.greedy;

import graph.Graph;

import java.util.Random;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: heuristics.greedy
 *         Description: The greedy degree heurisitc. Selects the vertex with
 *         minimum degree in the given graph.
 */
public class GreedyDegree implements GreedyCriteria {
  private Random rnd;
  private long seed;

  // This greedy class selects vertices by increasing degree.
  private final String criteria = "DEGREE";

  public GreedyDegree() {
    this(new Random().nextLong());
  }

  public GreedyDegree(long seed) {
    this.seed = seed;
    this.rnd = new Random(seed);
  }

  public GreedyDegree(Random random, long seed) {
    this.seed = seed;
    this.rnd = random;
  }

  public long getSeed() {
    return this.seed;
  }

  public Random getRnd() {
    return this.rnd;
  }

  public int getNextVertex(Graph G) {
    int minDegree = -1;
    int minVertex = -1;
    for (int i = 0; i < G.size(); i++) {
      int vertex = G.permutation[i];
      if (!G.contains(vertex))
        continue;

      if (G.degree(vertex) < minDegree || minDegree < 0) {
        minDegree = G.degree(vertex);
        minVertex = vertex;
      }
    }
    return minVertex;
  }

  public String getCriteriaAsString() {
    return this.criteria;
  }
}
