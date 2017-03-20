package heuristics.greedy;

import graph.Graph;

import java.util.ArrayList;
import java.util.Random;
import java.util.Set;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: heuristics.greedy
 *         Description: The min-fill heurisitc. Selects the vertex with
 *         least amount of fill-in.
 */
public class GreedyFillIn implements GreedyCriteria {
  private Random rnd;
  private long seed;

  // This greedy class selects vertices by first calculating the fill-in for
  // each vertex, and then selecting the vertex with the least amount of
  // fill-in required.
  private final String criteria = "FILLIN";

  public GreedyFillIn() {
    this(new Random().nextLong());
  }

  public GreedyFillIn(long seed) {
    this.seed = seed;
    this.rnd = new Random(seed);
  }

  public GreedyFillIn(Random random, long seed) {
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
    int minFillIn = -1;
    int minFillVertex = -1;
    for (int i = 0; i < G.size(); i++) {
      int vertex = G.permutation[i];
      if (!G.contains(vertex))
        continue;

      ArrayList<Integer> neighbours = new ArrayList<Integer>(G.neighboursOf(vertex));
      int numNeighbours = neighbours.size();
      int totalDeg = 0;

      int x, y;
      for (int j = 0; j < numNeighbours; j++) {
        x = neighbours.get(j);
        for (int k = 0; k < j; k++) {
          y = neighbours.get(k);
          if (G.isNeighbour(x, y))
            totalDeg++;
        }
      }
      int fillIn = ((numNeighbours*(numNeighbours - 1))/2) - totalDeg;

      // We will not find anything better.
      if (fillIn == 0) {
        return vertex;
      }

      if (fillIn < minFillIn || minFillIn < 0) {
        minFillIn = fillIn;
        minFillVertex = vertex;
      }
    }

    return minFillVertex;
  }

  public String getCriteriaAsString() {
    return this.criteria;
  }
}
