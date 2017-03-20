package graph;

import java.util.*;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: graph
 *         Description: A graph containing a set of vertices from 0..n-1.
 *         Stores the graph in both an adjacency list and adjacency matrix
 *         representation.
 */
public class Graph {
  // A boolean array, indexed from 0..n-1. Index i is true if vertex i is in
  // the graph.
  private boolean[] vertices;
  private HashSet<Integer>[] edges;
  private boolean[][] adjacencyMat;
  private String[] names;
  private int n;
  private int currentNumOfVertices;

  // The fill-in for every vertex is stored here. Used by the turbo-fillin
  // algorithm to sort by fill-in.
  public int[] fillIn;

  // Create a random permutation on the numbers 0..n-1. Used by the greedy
  // algorithms to break ties.
  public int[] permutation;

  public Graph(int n) {
    vertices = new boolean[n];
    edges = new HashSet[n];
    names = new String[n];
    fillIn = new int[n];
    this.n = n;
    this.currentNumOfVertices = 0;

    this.adjacencyMat = new boolean[n][n];

    permutation = new int[n];
    for (int i = 0; i < n; i++) {
      permutation[i] = i;
      names[i] = Integer.toString(i);
    }
  }

  /**
   * Generates a random permutation using the Fisher-Yates shuffle in O(n) time.
   * @param rnd The random object.
   */
  public void generateRandomPermutation(Random rnd) {
    for (int i = permutation.length - 1; i > 0; i--)
    {
      int index = rnd.nextInt(i + 1);
      int temp = permutation[index];
      permutation[index] = permutation[i];
      permutation[i] = temp;
    }
  }

  public void addEdge(int u, int v) {
    if (!contains(u) || !contains(v))
      return;

    edges[u].add(v);
    edges[v].add(u);
    this.adjacencyMat[u][v] = true;
    this.adjacencyMat[v][u] = true;
  }

  public void removeEdge(int u, int v) {
    if (!contains(u) || !contains(v))
      return;

    edges[u].remove(v);
    edges[v].remove(u);
    this.adjacencyMat[u][v] = false;
    this.adjacencyMat[v][u] = false;
  }

  // Removes the neighbour v from u.
  public void removeNeighbour(int u, int v) {
    if (!contains(u) || !contains(v))
      return;

    edges[u].remove(v);
    this.adjacencyMat[u][v] = false;
  }

  public void addVertex(int v) {
    vertices[v] = true;
    currentNumOfVertices++;

    if (edges[v] == null) {
      edges[v] = new HashSet<Integer>();
    }
  }

  public void addVertex(int v, String name) {
    vertices[v] = true;
    names[v] = name;
    currentNumOfVertices++;

    if (edges[v] == null) {
      edges[v] = new HashSet<Integer>();
    }
  }

  public void removeVertex(int v) {
    if (!contains(v))
      return;

    for (int neighbour : neighboursOf(v)) {
      edges[neighbour].remove(v);
    }

    // Setting this to false means we'll never read any information about it.
    vertices[v] = false;
    currentNumOfVertices--;
  }

  public Set<Integer> neighboursOf(int v) {
    if (!contains(v))
      return null;

    return edges[v];
  }

  public boolean isNeighbour(int v, int u) {
    if (!contains(u) || !contains(v))
      return false;

    return adjacencyMat[u][v] && adjacencyMat[v][u];
  }

  public int degree(int v) {
    if (!contains(v))
      return 0;
    return edges[v].size();
  }

  public int size() {
    return n;
  }

  // Inefficient method to determine the number of edges in the graph. Only
  // ever called for debugging purposes.
  public int numEdges() {
    int total = 0;
    for (int i = 0; i < this.size(); i++) {
      if (!this.contains(i))
        continue;
      total += this.neighboursOf(i).size();
    }
    return total/2;
  }

  public int getCurrentNumOfVertices() {
    return currentNumOfVertices;
  }

  public boolean contains(int v) {
    return vertices[v];
  }

  public String getVertexName(int v) {
    if (!vertices[v])
      return null;
    return names[v];
  }

  public boolean[] getVertices() {
    return vertices;
  }

  public int getVertexForName(String name) {
    for (int i = 0; i < names.length; i++) {
      if (name.equals(names[i]))
        return i;
    }
    return -1;
  }

  public String toString() {
    String out = "{";
    for (int i = 0; i < n; i++) {
      if (!contains(i))
        continue;
      String vertexString = getVertexName(i) + ": {";
      for (int j : neighboursOf(i)) {
        vertexString += getVertexName(j) + ", ";
      }
      if (neighboursOf(i).size() > 0)
        vertexString = vertexString.substring(0, vertexString.length()-2);
      vertexString += "}";
      out += vertexString + ", ";
    }
    if (size() > 0)
      out = out.substring(0, out.length()-2);
    return out + "}";
  }

  // Deep copies a graph in O(n + m) time.
  public Graph clone() {
    int n = size();
    Graph H = new Graph(n);
    for (int i = 0; i < n; i++) {
      H.permutation[i] = this.permutation[i];
      if (!contains(i))
        continue;
      H.addVertex(i, getVertexName(i));
    }

    for (int i = 0; i < n; i++) {
      if (!contains(i))
        continue;
      for (int j : neighboursOf(i))
        H.addEdge(i, j);
    }
    return H;
  }

  public static Graph computeFillInGraph(Graph G, ArrayList<Integer> ordering) {
    // Clone graph in O(n + m).
    Graph fillIn = G.clone();

    // For each vertex in ordering, go through all pairs of neighbours and
    // add edges.
    for (int i = 0; i < ordering.size(); i++) {
      int vertex = ordering.get(i);
      ArrayList<Integer> neighbours =
          new ArrayList<Integer>(fillIn.neighboursOf(vertex));
      for (int x : neighbours) {
        for (int y : neighbours) {
          if (x == y
              || ordering.indexOf(x) < i
              || ordering.indexOf(y) < i) {
            continue;
          }

          if (!fillIn.isNeighbour(x, y))
            fillIn.addEdge(x, y);
        }
      }
    }

    return fillIn;
  }

  public static void computeFillIn(Graph G, ArrayList<Integer> W) {
    for (int v : W) {
      ArrayList<Integer> neighbours = new ArrayList<Integer>(G.neighboursOf(v));
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
      G.fillIn[v] = fillIn;
    }
  }

  /**
   * Sorts a list of vertices W by a given criteria. Either "degree" or
   * "fillin". Degree will sort vertices by degree, from smallest to
   * largest. Fill-in will sort vertices by fill-in, with the vertex with the
   * least amount of fill-in appearing first.
   * @param W The set W to sort.
   */
  public void sortByCriteria(ArrayList<Integer> W, String criteria) {
    class VertexAttr implements Comparable<VertexAttr> {
      public Graph G;
      public int vertex;
      public String criteria;

      private String DEGREE_STR = "DEGREE";
      private String FILLIN_STR = "FILLIN";

      public VertexAttr(Graph G, int vertex, String criteria) {
        this.G = G;
        this.vertex = vertex;
        this.criteria = criteria;
      }

      @Override
      public int compareTo(VertexAttr other) {
        if (this.criteria.equalsIgnoreCase(DEGREE_STR)) {
          return Integer.compare(G.degree(this.vertex),
                                 G.degree(other.vertex));
        } else if (this.criteria.equalsIgnoreCase(FILLIN_STR)) {
          return Integer.compare(G.fillIn[this.vertex],
                                 G.fillIn[other.vertex]);
        }
        return 0;
      }
    }

    ArrayList<VertexAttr> sorted = new ArrayList<VertexAttr>();
    for (int u : W) {
      sorted.add(new VertexAttr(this, u, criteria));
    }

    Collections.sort(sorted);
    W.clear();

    for (VertexAttr u : sorted) {
      W.add(u.vertex);
    }
  }
}
