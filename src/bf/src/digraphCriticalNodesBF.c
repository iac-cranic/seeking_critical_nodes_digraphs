//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#include "graph.h"
#include <getopt.h>
#include <string.h> 

#define USAGE "\nUsage: digraphCriticalNodesBF --graph <FILE_NAME> --format <TYPE> --stop <NODES_NUMBER>\n\n"\
			  "\tStarting from size 1, we remove incremental sets of critical nodes from the graph, we stop when:\n"\
			  "\t1) we reach a 0 connectivity value or 2) we reach the size defined by the <stop> parameter.\n"\
			  "\tThe procedure computes the critical nodes of the graph based on the formula:\n\n"\
			  "\t\tf(G) = \\sum_{i=1}^{l} \binom{C_i}{2}\n\n"\
			  "\tWhere C_i are the strongly connected components of G, and a node x is a critical node\n"\
			  "\tof G if it minimizes f(G\\x), i.e x = arg min_{v \\in V} f(G\\v) \n\n"\
			  "\t-s, --stop <NODES_NUMBER>\n"\
			  "\t\tMaximum size of the critical nodes set to remove, 0 is interpreted as all.\n"\
			  "\t-f, --graph <FILE_NAME>\n"\
			  "\t\tRead the graph from file <FILE_NAME>.\n"\
			  "\t-t, --format <TYPE>\n"\
			  "\t\tThe format of the file containing the graph, admitted types are:\n"\
			  "\t\t\t- 0 (Edges List) The first line of the file is expected to contain the number of nodes\n"\
			  "\t\t\t    and edges of the graph in the following format: \"# Nodes: <NUMBER> Edges: <NUMBER>\"  \n\n"


int main(int argc, char **argv){
	
    int file_type;						/* Format of the file containing the graph */
    char ch;                    		/* Option read by getopt_long */ 
	char *file_name;      				/* Name of the file containing the graph */   
	digraph_type *g;
	vertex_type max_k;
	vertex_type k;
	
	file_name = NULL;
	ch = 0;
	file_type = -1;
	k = 1;
	max_k = 0;
	g = NULL; 
		
	static struct option long_options[] ={
	    {"graph", required_argument, NULL, 'f'},
	    {"format", required_argument, NULL, 't'},
	    {"stop", required_argument, NULL, 's'},
		{"help", no_argument, NULL, 'h'},
	    {NULL, 0, NULL, 0}
	};
	
	/* getopt_long() returns the option character when a short option is recognized. 
	 * For a long option, it returns <val> if <flag> is NULL, and 0 otherwise. */
	while ((ch = getopt_long(argc, argv, "f:t:s:h", long_options, NULL)) != -1){
	    switch (ch){
		    case 'f':
	            file_name = strdup(optarg);
	            break;
		    case 's':
	            max_k = atoi(optarg);
	            break;
			case 't':
	            file_type = atoi(optarg); /* Check that the file format is one of the allowed types */
				if ( file_type != EDGE_LIST ){
	                printf(USAGE);
	                exit(EXIT_FAILURE);
				}
	            break;
			case 'h':
            default:
                printf(USAGE);
                exit(EXIT_FAILURE);
                break;
	    }
	}

	
	/* Check that all needed information have been provided through the command line */
	if ( file_type == -1 || file_name == NULL) {
        printf(USAGE);
        exit(EXIT_FAILURE);
	}
	
	/* Read the graph */
	g = get_digraph_from_file( file_name, EDGE_LIST );
	if (g == NULL){ exit(EXIT_FAILURE); }
    free(file_name);	

	/* Remove Critical Nodes */
	uint64_t end_conn, init_conn, tot_sets;
	vertex_type **critical_sets;
	
	if (max_k == 0){
		max_k = g->nodes;
	}
	
    cc_type *scc;
    scc = get_digraph_scc(g);
    init_conn = get_pairs_connectivity(scc);
	fprintf(stdout,"%"PRIvertex ",%"PRIu64 "\n", 0, init_conn);
    free_cc(scc);
    end_conn = init_conn;
	while( end_conn != 0 && k <= max_k ){
		critical_sets = get_digraph_critical_nodes_bruteforce( g, k, 0, &end_conn, &tot_sets);
        free(critical_sets);
		fprintf(stdout,"%"PRIvertex ",%"PRIu64 "\n", k, end_conn);
		k++;
	}

	free_digraph(g);
}












