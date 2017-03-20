package graph;

import java.util.ArrayList;

/**
 * @author mitchelljones
 *         Project: tree_decomp
 *         Package: graph
 *         Description: The vertex/bag in a tree decomposition.
 */
public class TreeVertex {
  private String name;
  private ArrayList<Integer> bag;
  private ArrayList<TreeVertex> neighbours;

  // Record which vertex was eliminated in the graph when constructed this
  // vertex in the tree decomposition.
  private int eliminatedVertex;

  public TreeVertex(String name, int eliminatedVertex, ArrayList<Integer> bag) {
    this.name = name;
    this.bag = bag;
    this.eliminatedVertex = eliminatedVertex;
    this.neighbours = new ArrayList<TreeVertex>();
  }

  public TreeVertex(String name, int eliminatedVertex) {
    this(name, eliminatedVertex, new ArrayList<Integer>());
  }

  /**
   * @return The name of the tree vertex.
   */
  public String getName() {
    return this.name;
  }

  /**
   * @return The vertex that was eliminated from the graph when constructing
   * this vertex in the tree decomposition.
   */
  public int getEliminatedVertex() {
    return this.eliminatedVertex;
  }

  /**
   * @return The size of the bag.
   */
  public int sizeOfBag() {
    return this.bag.size();
  }

  /**
   * Adds a vertex u to the bag.
   * @param u
   */
  public void addToBag(int u) {
    this.bag.add(u);
  }

  /**
   * @return The bag of vertices in the tree vertex.
   */
  public ArrayList<Integer> bag() {
    return this.bag;
  }

  /**
   * Adds a neighbour to the tree vertex.
   * @param u The new neighbour.
   */
  public void addNeighbour(TreeVertex u) {
    this.neighbours.add(u);
  }

  /**
   * @return A collection of the neighbouring tree vertices in the tree
   * decomposition.
   */
  public ArrayList<TreeVertex> getNeighbours() {
    return this.neighbours;
  }

  /**
   * Removes a neighbour u from this vertex.
   * @param u The neighbour to be removed.
   */
  public void removeNeighbour(TreeVertex u) {
    this.neighbours.remove(u);
  }

  /**
   * @param u The neighbour to test.
   * @return True if u is a neighbour of this vertex, false otherwise.
   */
  public boolean isNeighbour(TreeVertex u) {
    return this.neighbours.contains(u);
  }

  public String toString() {
    String name = "(";
    for (int v : this.bag) {
      name += v + ",";
    }
    if (this.bag.size() > 0)
      name = name.substring(0, name.length()-1);
    name += ")";

    String out = name + ": {";
    for (TreeVertex u : this.neighbours) {
      String neighbourName = "(";
      for (int v : u.bag()) {
        neighbourName += v + ",";
      }
      if (this.bag.size() > 0)
        neighbourName = neighbourName.substring(0, neighbourName.length()-1);
      neighbourName += ")";
      out += neighbourName + ", ";
    }
    if (this.neighbours.size() > 0)
      out = out.substring(0, out.length()-2);
    return out + "}";
  }
}
