package heuristics.greedy;

import graph.Graph;

import java.util.Random;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: heuristics.greedy
 *         Description: An interface all greedy algorithms must implement.
 *         Allows for the customisation of selecting the next node to
 *         add to the ordering.
 */
public interface GreedyCriteria {
  int getNextVertex(Graph G);
  long getSeed();
  Random getRnd();

  // Each class that implements GreedyCriteria must return a unique string,
  // which describes the greedy heuristic. For example, "DEGREE" denotes the
  // greedy method that selects vertices with the smallest degree first. This
  // is used by the Turbocharge class, when sorting vertices by the greedy
  // criteria. See method Graph.sortByCriteria.
  String getCriteriaAsString();
}
