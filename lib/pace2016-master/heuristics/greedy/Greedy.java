package heuristics.greedy;

import graph.Graph;
import heuristics.ConstructTree;

import java.util.ArrayList;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: heuristics.greedy
 *         Description: Runs a greedy algorithm given the criteria and graph.
 */
public class Greedy {

  /**
   * Runs the GreedyX algorithm. Using the given a criteria to select the
   * next node, it constructs the fill-in graph H and returns an elimination
   * ordering. This method breaks ties using the `criteria` Random object.
   * @param G The graph G.
   * @param criteria The criteria.
   * @return An ordering of vertices, represented by their index.
   */
  public static ArrayList<Integer> greedyWithCriteria(Graph G,
                                                     GreedyCriteria criteria) {
    G.generateRandomPermutation(criteria.getRnd());
    Graph H = G.clone();
    ArrayList<Integer> ordering = new ArrayList<Integer>();
    while (H.getCurrentNumOfVertices() > 0) {
      int v = criteria.getNextVertex(H);
      ordering.add(v);
      ConstructTree.eliminateVertex(H, v);
    }
    return ordering;
  }
}
