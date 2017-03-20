#pragma once
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <vector>

#include "always_assert.hpp"
#include "graph.hpp"

class TD {
 public:
  using Bag = size_t;

  TD() {}
  TD(std::vector<Vertex> &&starting_bag) {
    width_ = starting_bag.size();
    bags.emplace_back(starting_bag);
    parent.push_back(0);
  }

  Bag add_leaf(Vertex v) {
    bags.push_back({v});
    parent.push_back(0);
    width_ = std::max(size_t(1), width_);
    return bags.size() - 1;
  }

  Bag add_child(Bag p, std::vector<Vertex> &&contents) {
    always_assert(0 <= p && p < bags.size());
    width_ = std::max(contents.size(), width_);
    bags.emplace_back(contents);
    parent.push_back(p);
    return bags.size() - 1;
  }

  size_t width() const { return width_; }

  std::string to_string(const Graph &g) const {
    std::stringstream ss;

    ss<< "s td " << bags.size() << ' ' << width() << ' '
        << g.num_vertices() << "\n";
    for (size_t i = 0; i < bags.size(); ++i) {
      ss << "b " << i + 1 << ' ';
      for (Vertex v : bags[i])
        ss << v + 1 << ' ';
      ss << "\n";
    }
    for (size_t i = 1; i < bags.size(); ++i)
      ss << parent[i] + 1 << ' ' << i + 1 << "\n";

    return ss.str();
  }

  void print(const Graph &g, std::ostream &out = std::cout) const {
    out << to_string(g);
  }

  void swap(TD &other) {
    bags.swap(other.bags);
    parent.swap(other.parent);
    std::swap(width_, other.width_);
  }

  bool is_valid(const Graph &graph) const {
    std::vector<Vertex> vertices_seen(graph.num_vertices(), false);

    std::set<Edge> all_edges;
    for (Vertex v : graph.vertices()) 
      for (Vertex w : graph.neighbors(v))
        all_edges.insert({v, w});

    for (const BagImpl &bag : bags) {
      for (Vertex v : bag) {
        vertices_seen[v] = true;
        for (Vertex w : bag) {
          all_edges.erase({v, w});
        }
      }
    }

    if (!all_edges.empty())
      return false;

    for (bool b : vertices_seen)
      if (!b)
        return false;

    if (parent.size() < 2)
      return true;

    always_assert(parent.size() == graph.num_vertices());
    std::vector<bool> forgotten(graph.num_vertices(), false);

    for (size_t i = graph.num_vertices() - 1; i > 0; --i) {
      if (parent[i] >= i)
        return false;

      std::set<Vertex> in_parent(bags[parent[i]].begin(), bags[parent[i]].end());

      for (Vertex v : bags[i]) {
        if (forgotten[v])
          return false;
        if (in_parent.count(v) == 0)
          forgotten[v] = true;
      }
    }

    return true;
  }

 private:
  using BagImpl = std::vector<Vertex>;
  std::vector<BagImpl> bags;
  std::vector<Bag> parent;
  size_t width_ = 0;
};
