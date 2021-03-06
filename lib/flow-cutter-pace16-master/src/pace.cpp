
#include "id_func.h"
#include "list_graph.h"
#include "multi_arc.h"
#include "sort_arc.h"
#include "chain.h"
#include "flow_cutter.h"
#include "greedy_order.h"

#include "node_flow_cutter.h"
#include "contraction_graph.h"
#include "cch_order.h"
#include "tree_decomposition.h"
#include "separator.h"
#include <iostream>
#include <limits>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#ifdef PARALLELIZE
#include <omp.h>
#include <atomic>
#endif

#include <sys/time.h>
#include <unistd.h>
using namespace std;

ArrayIDIDFunc tail, head;
const char*volatile best_decomposition = 0;
int best_bag_size = numeric_limits<int>::max();

void ignore_return_value(int){}

int compute_max_bag_size(const ArrayIDIDFunc&order){
	auto inv_order = inverse_permutation(order);
	int current_tail = -1;
	int current_tail_up_deg = 0;
	int max_up_deg = 0;
	compute_chordal_supergraph(
		chain(tail, inv_order), chain(head, inv_order), 
		[&](int x, int y){
			if(current_tail != x){
				current_tail = x;
				max_to(max_up_deg, current_tail_up_deg);
				current_tail_up_deg = 0;
			}
			++current_tail_up_deg;
		}
	);
	return max_up_deg+1;
}

unsigned long long get_milli_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (unsigned long long)(tv.tv_sec) * 1000
	   + (unsigned long long)(tv.tv_usec) / 1000;
}

const char*compute_decomposition(const ArrayIDIDFunc&order){
	ostringstream out;
	print_tree_decompostion(out, tail, head, move(order));
	char*buf = new char[out.str().length()+1];
	memcpy(buf, out.str().c_str(), out.str().length()+1);
	return buf;
}

void test_new_order(ArrayIDIDFunc order){
	int x = compute_max_bag_size(order);
	#ifdef PARALLELIZE
	#pragma omp critical
	#endif
	{
		if(x < best_bag_size){
			best_bag_size = x;
			const char*old_decomposition = best_decomposition;
			best_decomposition = compute_decomposition(move(order));
			delete[]old_decomposition;
			{
				string msg = "c status "+to_string(best_bag_size)+" "+to_string(get_milli_time())+"\n";
				ignore_return_value(write(STDOUT_FILENO, msg.data(), msg.length()));
			}
		}
	}
}

char no_decomposition_message[] = "c info programm was aborted before any decomposition was computed\n";

#ifdef PARALLELIZE
volatile atomic_flag only_one_thread_in_signal_handler = ATOMIC_FLAG_INIT;
#endif

void signal_handler(int)
{
	#ifdef PARALLELIZE
	while (only_one_thread_in_signal_handler.test_and_set()) {}
	#endif

	const char*x = best_decomposition;
	if(x != 0)
		ignore_return_value(write(STDOUT_FILENO, x, strlen(x)));
	else
		ignore_return_value(write(STDOUT_FILENO, no_decomposition_message, sizeof(no_decomposition_message)));

	_Exit(EXIT_SUCCESS);
}

/*
int main(int argc, char * argv[]){
	string file_name = "-";
	if(argc == 2)
		file_name = argv[1];
	auto g = uncached_load_pace_graph(file_name);
	tail = std::move(g.tail);
	head = std::move(g.head);
	
	int random_seed = 0;
	
	if(argc == 3){
		if(string(argv[1]) == "-s"){
			random_seed = atoi(argv[2]);
		}
	}
	std::minstd_rand rand_gen;
	rand_gen.seed(
		random_seed 
	);

	flow_cutter::Config config;
	config.cutter_count = 1;
	config.random_seed = rand_gen();
	
	++config.cutter_count;
	test_new_order(cch_order::compute_cch_graph_order(tail, head, flow_cutter::ComputeSeparator(config)));
	*/
	/*
	for(int i=0;i<2;++i){
		config.random_seed = rand_gen();
		if(i % 32 == 0)
			++config.cutter_count;
		test_new_order(cch_order::compute_cch_graph_order(tail, head, flow_cutter::ComputeSeparator(config)));
		cout << "step3" << endl;
	}
	*/
	/*
	const char*x = best_decomposition;
	if(x != 0)
		ignore_return_value(write(STDOUT_FILENO, x, strlen(x)));
	else
		ignore_return_value(write(STDOUT_FILENO, no_decomposition_message, sizeof(no_decomposition_message)));
}
*/


int main(int argc, char*argv[]){
	signal(SIGTERM, signal_handler);
	signal(SIGINT, signal_handler);

	signal(SIGSEGV, signal_handler);
	try{
		{
			string file_name = "-";
			if(argc == 2)
				file_name = argv[1];
			auto g = uncached_load_pace_graph(file_name);
			tail = std::move(g.tail);
			head = std::move(g.head);
		}

		int random_seed = 0;

		if(argc == 3){
			if(string(argv[1]) == "-s"){
				random_seed = atoi(argv[2]);
			}
		}

		#ifdef PARALLELIZE
		#pragma omp parallel
		#endif
		{
			try{
				#ifdef PARALLELIZE
				#pragma omp sections nowait
				#endif
				
				{
					test_new_order(compute_greedy_min_degree_order(tail, head));
					#ifdef PARALLELIZE
					#pragma omp section
					#endif
					test_new_order(compute_greedy_min_shortcut_order(tail, head));
				}
				
				std::minstd_rand rand_gen;
				rand_gen.seed(
					random_seed 
					#ifdef PARALLELIZE
					+ omp_get_thread_num()
					#endif
				);

				flow_cutter::Config config;
				config.cutter_count = 1;
				config.random_seed = rand_gen();

				for(int i=0;;++i){
					config.random_seed = rand_gen();
					if(i % 32 == 0)
						++config.cutter_count;
			
					test_new_order(cch_order::compute_cch_graph_order(tail, head, flow_cutter::ComputeSeparator(config)));

				}
			}catch(...){
				signal_handler(0);
			}
		}
	}catch(...){
		signal_handler(0);
	}
}

