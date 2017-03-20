package heuristics;

import graph.*;

import java.util.ArrayList;
import java.util.Set;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: heuristics
 *         Description: Constructs a tree decomposition given an elimination
 *         ordering (permutation) on the vertices. Also contains utility
 *         methods such as eliminating a vertex from the graph.
 */
public class ConstructTree {

  /**
   * Constructs a tree decomposition given an elimination ordering
   * (permutation) on the nodes.
   * @param G The graph G.
   * @param ordering An ordering of the nodes to eliminate, given by their
   *                 unique name.
   * @return A tree decomposition.
   */
  public static TreeGraph orderingToTreeGraph(Graph G,
                                              ArrayList<Integer> ordering) {
    // Step 1. Create all of the vertices in the tree decomposition by
    // eliminating vertices in the graph.
    Graph H = G.clone();
    TreeGraph tree = new TreeGraph();
    for (int i = 0; i < ordering.size(); i++) {
      // Select the next vertex in the ordering, create a tree node and bag
      // (corresponding to the neighbours), and eliminate the vertex.
      int v = ordering.get(i);
      ArrayList<Integer> bag = new ArrayList<Integer>();
      for (int u : H.neighboursOf(v)) {
        bag.add(u);
      }
      bag.add(v);
      TreeVertex newNode = new TreeVertex(bag.toString(), v, bag);
      tree.addTreeVertex(newNode);
      eliminateVertex(H, v);
    }

    // Step 2. Compute fill-in graph given the ordering.
    Graph fillIn = Graph.computeFillInGraph(G, ordering);

    // Step 3. Add the edges between nodes in the tree decomposition.
    for (int i = 0; i < ordering.size(); i++) {
      // For each vertex v in the ordering, find the first neighbour u in
      // the fill in graph that appears after v. Connect the two tree nodes
      // that were created by eliminating v and u.
      int v = ordering.get(i);
      int lowestNeighbour = -1;
      int lowestNumNeighbour = ordering.size();
      for (int u : fillIn.neighboursOf(v)) {
        int index = ordering.indexOf(u);
        if (index < lowestNumNeighbour && index > i) {
          lowestNumNeighbour = index;
          lowestNeighbour = u;
        }
      }
      if (lowestNeighbour >= 0) {
        TreeVertex parent = tree.treeVertexCreatedByEliminating(v);
        TreeVertex child = tree.treeVertexCreatedByEliminating(lowestNeighbour);
        tree.addEdge(parent, child);
      }
    }

    return tree;
  }

  /**
   * Eliminates a vertex v from the graph G, and completes the neighbouring
   * subgraph to form a clique.
   * @param G The graph G.
   * @param vertex The vertex v to eliminate.
   * @return
   */
  public static void eliminateVertex(Graph G, int vertex) {
    ArrayList<Integer> neighbours =
        new ArrayList<Integer>(G.neighboursOf(vertex));
    int numNeighbours = neighbours.size();
    int x, y;
    for (int i = 0; i < numNeighbours; i++) {
      x = neighbours.get(i);
      for (int j = 0; j < i; j++) {
        y = neighbours.get(j);
        if (!G.isNeighbour(x, y))
          G.addEdge(x, y);
      }
    }
    G.removeVertex(vertex);
  }
}
