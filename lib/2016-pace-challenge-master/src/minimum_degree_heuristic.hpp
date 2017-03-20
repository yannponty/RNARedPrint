#pragma once
#include "tree_decomposition.hpp"
#include "graph.hpp"

TD minimum_degree_heuristic(Graph graph, size_t lb = std::numeric_limits<size_t>::max());
TD minimum_fillin_heuristic(Graph graph, size_t lb = std::numeric_limits<size_t>::max());
