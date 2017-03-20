#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "always_assert.hpp"
#include "graph.hpp"
#include "minimum_degree_heuristic.hpp"
#include "tree_decomposition.hpp"

namespace {
std::string *tmp_str(new std::string(""));
std::atomic<std::string *> best_td_str(new std::string(""));
std::atomic<size_t> best_width(0);

void signal_handler(int signum) {
  if (signum == SIGTERM) {
    std::cout << *best_td_str.load();
    std::exit(0);
  }
  std::cerr << "Invalid signal, aborting\n";
  std::exit(4);
}

#ifdef VALIDATE_TD
void validate_td(const Graph &graph, const TD &td) {
  always_assert(td.is_valid(graph));
}
#else
void validate_td(const Graph &, const TD &) {}
#endif

}  // anonymous namespace

int main(int argc, char **argv) {
  // Set up signal handling
  struct sigaction sa;
  sa.sa_handler = signal_handler;
  sa.sa_flags = SA_RESTART;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGTERM);
  sigaction(SIGTERM, &sa, NULL);

  int opt;
  while ((opt = getopt(argc, argv, "fds:")) != -1) {
    if (opt == 's') {
      srand(unsigned(std::stoul(optarg)));
    } else {
      std::cerr << "Invalid argument: " << opt << ", aborting\n";
      return 1;
    }
  }

  Graph graph;

  if (optind == argc - 1) {
    std::ifstream graph_input(argv[argc - 1]);
    graph = Graph(graph_input);
  } else {
    graph = Graph(std::cin);
  }

  TD td = minimum_degree_heuristic(graph);
  validate_td(graph, td);

  *tmp_str = td.to_string(graph);
  tmp_str = best_td_str.exchange(tmp_str);
  best_width = td.width();
  std::cout << "c status " << td.width() << " "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::system_clock::now().time_since_epoch()).count()
            << "\n";
            
	td = minimum_fillin_heuristic(graph, best_width.load());
  validate_td(graph, td);
	
  while (true) {
    td = minimum_fillin_heuristic(graph, best_width.load());
    validate_td(graph, td);
    if (td.width() < best_width.load()) {
      *tmp_str = td.to_string(graph);
      tmp_str = best_td_str.exchange(tmp_str);
      best_width = td.width();
      std::cout << "c status " << td.width() << " "
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       std::chrono::system_clock::now().time_since_epoch())
                       .count() << "\n"; 
    }
  }
 
}
