//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <limits.h>
#include <igraph.h>
#include <signal.h>
#include "define.h"
#include "wrapper_igraph.h"
#include "timer.h"



#define USAGE "\nUsage: igraph_random --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME>\n\n"\
	"\tThe program tests the auxiliary wrapper functions to igraph\n\n"\
"\t-f, --graph <FILE_NAME>\n"\
"\t\tRead the graph from file <FILE_NAME>.\n"\
"\t-o, --outfile <FILE_NAME>\n"\
"\t\tThe output filename <FILE_NAME>.\n"\
"\t-p, percentage <PERCENTAGE>.\n"\
"\t\tThe percentage of nodes to output.\n"\
"\t-t, --format <TYPE>\n"\
"\t\tThe format of the file containing the graph, admitted types are:\n"\
"\t\t\t- 0 (DIMACS)\n"\
"\t\t\t- 1 (Edges List) The first line of the file is expected to contain the number of nodes\n"\
"\t\t\t    and edges of the graph in the following format: \"# Nodes: <NUMBER> Edges: <NUMBER>\"  \n\n"

FILE *fout=NULL;

// Type of graph files supported
const char *GTYPE_NAMES[2] = {"DIMACS", "RMAT"};

void sigUser1Handler(int sig){
	if(NULL != fout){
		fclose(fout);
	}
	fprintf(stderr,"[WARN] SIGUSR1 received - closing\n");
	exit(EXIT_FAILURE);
}


int main(int argc, char ** argv){

    igraph_t graph;
    igraph_vector_t v;
    igraph_vector_t result;
    VTYPE root;
    igraph_real_t *edges;


    struct sigaction new_act;
    int         gtype;                  // Graph type (see enum graph_type)
    char        c;                      // getopt char
    int 	    read_edges;				// Determine if read_graph should read also the edges list
    VTYPE       n, m;                  // n -> number of vertices, m -> number of edges
    VTYPE	    node_perc;		       // percentage of node to output
    // Others
    char *graph_file_name;             // input file name
    char *out_file_name;
    int i ;
    int metric;

    // Initialize
    gtype     = -1;
    read_edges=  0;
    n         =  0;
    m         =  0;
    node_perc =  0;

    graph_file_name = NULL;
    out_file_name = NULL;
    i = 0;
    metric = 0;
    fout = NULL;


    static struct option long_options[] ={
        {"graph", required_argument, NULL, 'f'},
        {"format", required_argument, NULL, 't'},
        {"outfile", required_argument, NULL, 'o'},
        {"percentage", required_argument, NULL, 'p'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "f:t:m:o:p:h", long_options, NULL)) != EOF) {
        switch (c) {
            // Graph file name
            case 'f':
                graph_file_name = StringCopy(optarg);
                break;
	    // Metric
	    case 'm':
		metric = atoi(optarg);
		break;
            // Output file name
            case 'o':
                out_file_name = StringCopy(optarg);
                break;
             // Percentage of node to output
            case 'p':
                node_perc = atoi(optarg);
                break;
            // Graph file type
            case 't':
                gtype = atoi(optarg);
                break;
            // Help
            case 'h':       
                fprintf(stderr,"%s",USAGE);
                exit(EXIT_FAILURE);
                break;
            default:
                fprintf(stderr,"%s",USAGE);
                exit(EXIT_FAILURE);
                break;
        }
    }

    if (gtype) {
        if (graph_file_name == NULL) {
            fprintf(stderr,"%s",USAGE);
            fprintf(stderr, "option -t require -f <input file>\n");
            exit(EXIT_FAILURE);
        }
    }
    if (graph_file_name) {
        if (gtype < 0) {
            fprintf(stderr,"%s",USAGE);
            fprintf(stderr, "option -f require -t <graph type>\n");
            exit(EXIT_FAILURE);
        }
    }
    if (gtype > RMAT) {
        fprintf(stderr,"%s",USAGE);
        fprintf(stderr, "file type not recognized!\n");
        exit(EXIT_FAILURE);
    }
    if( out_file_name == NULL){
        fprintf(stderr,"%s",USAGE);
        fprintf(stderr,"option -o <output filename> is mandatory\n");
        exit(EXIT_FAILURE);
    }
    if( node_perc <= 0 || node_perc > 100){
        fprintf(stderr,"%s",USAGE);
        fprintf(stderr,"option -p is mandatory and it should be 0 < p <= 100\n");
        exit(EXIT_FAILURE);
    }
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    edges = read_graph_igraph(gtype, graph_file_name, &n, &m, &root, 1);
    Free(graph_file_name);
    igraph_vector_view(&v, edges, 2*m);
    igraph_create(&graph, &v, 0, IGRAPH_DIRECTED);

    
    char buf[BUFFSIZE];
    for( i = 0; i < n ; i++){
        SETVAN(&graph, ID_STRING, i, i+1);
    }

    igraph_vector_init(&result,0);

    igraph_vector_t id_vector;

    igraph_vector_init(&id_vector, igraph_vcount(&graph));


    unsigned long long con = getConnectivity(&graph);

    // Set up sighandler for SIGUSR1
    new_act.sa_handler=sigUser1Handler;
    sigemptyset(&new_act.sa_mask);
    new_act.sa_flags=0;
    sigaction(SIGUSR1,&new_act,NULL);
    int ret = 0;

    unsigned int isolated = 0;
    unsigned int * removed = Calloc(graph.n, sizeof(int));
    unsigned int removed_num = 0;
    unsigned int max_node_to_remove = (n / 100.0) * node_perc;

    // Generate a pseudo-random permutation to be used for randomly delete vertices
    srand48(time(NULL));

    fout = Fopen(out_file_name, "w+");    
    Free(out_file_name);
    TIMER_DEF;
    TIMER_START;
    fprintf(fout,"%f,%d,%lu\n",0.0,removed_num,con);
    while ( con > 0 && removed_num < max_node_to_remove){

	VTYPE id_best=lrand48() % graph.n;
        TIMER_STOP;
	fprintf(fout,"%f,%d,",(TIMER_ELAPSED /  1000000.0), (int)VAN(&graph, ID_STRING, id_best));
	removed[removed_num] = VAN(&graph, ID_STRING, id_best);
        removed_num++;	

        ret = igraph_delete_vertices(&graph,igraph_vss_1(id_best));
        if ( ret == IGRAPH_EINVVID){
            fprintf(stderr,"[ERROR]: unable to remove %d\n",VAN(&graph, ID_STRING,id_best));
        }	

	con = getConnectivity(&graph);
	fprintf(fout,"%lu\n",con);

    }

    igraph_destroy(&graph);
    Free(removed);
    fclose(fout);
    return 0;
}


