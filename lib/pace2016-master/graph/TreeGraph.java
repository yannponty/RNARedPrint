package graph;

import java.util.*;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: graph
 *         Description: A graph representing a tree decomposition.
 */
public class TreeGraph {
  private ArrayList<TreeVertex> vertices;

  // Stores a hashmap of Integer (Vertex) -> TreeVertex pairs. Given an
  // integer representing a vertex v, it returns the TreeVertex that was
  // constructed by eliminating v.
  private HashMap<Integer, TreeVertex> eliminations;

  // The tree width plus one.
  private int treeWidthPlusOne;

  private int edges;

  public TreeGraph() {
    this.vertices = new ArrayList<TreeVertex>();
    this.eliminations = new HashMap<Integer, TreeVertex>();
    treeWidthPlusOne = 0;
    edges = 0;
  }

  /**
   * Adds an edge (v,u) to the tree decomposition.
   * @param v The first tree vertex.
   * @param u The second tree vertex.
   */
  public void addEdge(TreeVertex v, TreeVertex u) {
    if (v.isNeighbour(u))
      return;

    v.addNeighbour(u);
    u.addNeighbour(v);
    edges++;
  }

  /**
   * Returns the number of vertices in the tree decomposition.
   * @return The number of vertices.
   */
  public int size() {
    return this.vertices.size();
  }

  public int edges() {
    return this.edges;
  }

  /**
   * @param t The tree vertex.
   * @return The bag of vertices for a given tree vertex t.
   */
  public ArrayList<Integer> bagForTreeVertex(TreeVertex t) {
    return t.bag();
  }

  /**
   * Adds a tree vertex to the graph.
   * @param t The new tree vertex.
   */
  public void addTreeVertex(TreeVertex t) {
    if (t.sizeOfBag() > this.treeWidthPlusOne)
      this.treeWidthPlusOne = t.sizeOfBag();
    this.vertices.add(t);
    this.eliminations.put(t.getEliminatedVertex(), t);
  }

  /**
   * @return A collection of all tree vertices in the graph.
   */
  public ArrayList<TreeVertex> getTreeVertices() {
    return this.vertices;
  }

  /**
   * Returns the tree vertex that was created by eliminating the given vertex
   * v. Used in the greedy heuristics where there is a specific elimination
   * ordering.
   * @param v The vertex v.
   * @return The tree vertex created by eliminating v.
   */
  public TreeVertex treeVertexCreatedByEliminating(int v) {
    return eliminations.get(v);
  }

  /**
   * @return The tree width of the tree decomposition.
   */
  public int width() {
    return this.treeWidthPlusOne - 1;
  }

  public String toString() {
    String out = "{";
    for (TreeVertex v : this.vertices) {
      out += v.toString() + ", ";
    }
    if (this.vertices.size() > 0)
      out = out.substring(0, out.length()-2);
    return out + "}";
  }

  /**
   * This method finds all connected components in the tree. Used to ensure
   * that the final output is a connected tree, and not a forest.
   * @return A set of the connect components in the tree.
   */
  public ArrayList<ArrayList<TreeVertex>> findComponents() {
    ArrayList<ArrayList<TreeVertex>> components =
        new ArrayList<ArrayList<TreeVertex>>();
    Set<TreeVertex> visited = new HashSet<TreeVertex>();
    ArrayList<TreeVertex> queue = new ArrayList<TreeVertex>();

    while (visited.size() < this.size()) {
      // Find a vertex that has not been visited, and explore it's components.
      ArrayList<TreeVertex> component = new ArrayList<TreeVertex>();
      for (TreeVertex v : this.getTreeVertices()) {
        if (!visited.contains(v)) {
          queue.add(v);
          break;
        }
      }

      // BFS to find all vertices in this component.
      while (!queue.isEmpty()) {
        TreeVertex v = queue.remove(0);
        visited.add(v);
        component.add(v);
        for (TreeVertex neighbour : v.getNeighbours()) {
          if (!visited.contains(neighbour))
            queue.add(neighbour);
        }
      }

      components.add(component);
    }

    return components;
  }

  /**
   * An O(n) algorithm that checks if the tree decomposition contains a
   * cycle, using a modified DFS algorithm. Can be used to assert that a
   * tree decomposition construction algorithm correctly constructs a tree.
   * @return true if the tree decomposition contains a cycle.
   */
  public boolean containsCycle() {
    HashMap<TreeVertex, Integer> visitTime = new HashMap<TreeVertex, Integer>();
    Stack<TreeVertex> stack = new Stack<TreeVertex>();
    HashMap<TreeVertex, Iterator<TreeVertex>> neighbours =
        new HashMap<TreeVertex, Iterator<TreeVertex>>();
    for (TreeVertex v : this.getTreeVertices()) {
      neighbours.put(v, v.getNeighbours().iterator());
    }
    int time = 0;

    while (visitTime.size() < this.size()) {
      // Find the next node that hasn't been visited if there are multiple
      // components.
      for (TreeVertex v : this.getTreeVertices()) {
        if (!visitTime.containsKey(v)) {
          stack.push(v);
          break;
        }
      }

      while (!stack.empty()) {
        TreeVertex v = stack.peek();
        if (!visitTime.containsKey(v)) {
          visitTime.put(v, time++);
        }
        // If we find a neighbour that has been visited previously, that is
        // not v's direct parent (i.e. next on top of the stack) then G
        // contains a cycle.
        if (neighbours.get(v).hasNext()) {
          TreeVertex u = neighbours.get(v).next();
          if (!visitTime.containsKey(u)) {
            stack.push(u);
          } else {
            v = stack.pop();
            if (u != stack.peek() && visitTime.get(u) < visitTime.get(v)) {
              return true;
            } else {
              stack.push(v);
            }
          }
        }
        else {
          stack.pop();
        }
      }
    }

    return false;
  }
}
