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



#define BETWEENNESS 1
#define CLOSENESS_ALL 2
#define CLOSENESS_IN 21
#define CLOSENESS_OUT 22
#define DEGREE_ALL 3
#define DEGREE_IN 31
#define DEGREE_OUT 32
#define PAGERANK 4

#define BETWEENNESS_CMD(graph,result) igraph_betweenness(&graph, &result, igraph_vss_all(), IGRAPH_DIRECTED,0, 1)
#define CLOSENESS_ALL_CMD(graph,result) igraph_closeness(&graph, &result, igraph_vss_all(), IGRAPH_ALL, 0, 1)
#define CLOSENESS_IN_CMD(graph,result) igraph_closeness(&graph, &result, igraph_vss_all(), IGRAPH_IN, 0, 1)
#define CLOSENESS_OUT_CMD(graph,result) igraph_closeness(&graph, &result, igraph_vss_all(), IGRAPH_OUT, 0, 1)
#define DEGREE_ALL_CMD(graph,result) igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL, 0)
#define DEGREE_IN_CMD(graph,result) igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_IN, 0)
#define DEGREE_OUT_CMD(graph,result) igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_OUT, 0)
#define PAGERANK_CMD(graph,result) igraph_pagerank(&graph, IGRAPH_PAGERANK_ALGO_PRPACK, &result, 0, igraph_vss_all(), 0, 0.85, 0, 0);

#define USAGE "\nUsage: igraph_standard --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME> --metric <METRIC_TYPE> --percentage <PERCENTAGE>\n\n"\
                   "\t-f, --graph <FILE_NAME>\n"\
                       "\t\tRead the graph from file <FILE_NAME>.\n"\
                   "\t-m, --metric <METRIC_TYPS>\n"\
                       "\t\tSupported metrics are:\n"\
                          "\t\t\t- 1 Betwenness centrality\n"\
                          "\t\t\t- 2 Closeness (all)\n"\
                          "\t\t\t- 21 Closeness (in)\n"\
                          "\t\t\t- 22 Closeness (out)\n"\
                          "\t\t\t- 3 Degree (all)\n"\
                          "\t\t\t- 31 Degree (in)\n"\
                          "\t\t\t- 32 Degree (out)\n"\
                          "\t\t\t- 4 Pagerank\n"\
                   "\t-o, --outfile <FILE_NAME>\n"\
                       "\t\tThe name of the output file.\n"\
                   "\t-p, --percentage <PERCENTAGE>.\n"\
                       "\t\tThe percentage of nodes to remove.\n"\
                   "\t-t, --format <TYPE>\n"\
                       "\t\tThe format of the file containing the graph, supported types are:\n"\
                          "\t\t\t- 0 (DIMACS)\n"\
                          "\t\t\t- 1 (Edges List) The first line of the file is expected to contain the number of nodes\n"\
                          "\t\t\t    and edges of the graph in the following format: \"# Nodes: <NUMBER> Edges: <NUMBER>\"  \n\n"

FILE *fout=NULL;
FILE *ftime=NULL;

// Type of graph files supported
const char *GTYPE_NAMES[2] = {"DIMACS", "RMAT"};

void sigUser1Handler(int sig){
	if(NULL != fout){
		fclose(fout);
	}
	if(NULL != ftime){
		fclose(ftime);
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
    int         c;                      // getopt char
    int 	    read_edges;				// Determine if read_graph should read also the edges list
    VTYPE       n, m;                  // n -> number of vertices, m -> number of edges
    VTYPE	    node_perc;		       // percentage of node to output
    // Others
    char *graph_file_name;             // input file name
    char *out_file_name;
    char *time_file_name;
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
    time_file_name= NULL;
    i = 0;
    metric = 0;


    static struct option long_options[] ={
        {"graph", required_argument, NULL, 'f'},
        {"format", required_argument, NULL, 't'},
	{"metric", required_argument, NULL, 'm'},
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

    if(!metric){
       fprintf(stderr,"%s",USAGE);
       fprintf(stderr,"Invalid metric selected\n");
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
    TIMER_DEF;
    TIMER_START;
    if(!metric){
	    fprintf(stderr,"%s",USAGE);
	    fprintf(stderr,"option -m is mandatory\n");
	    exit(EXIT_FAILURE);
    } else {
	    switch(metric){
		    case BETWEENNESS:
			    ret=BETWEENNESS_CMD(graph,result);
			    break;
		    case CLOSENESS_ALL:
			    ret=CLOSENESS_ALL_CMD(graph,result);
			    break;
		    case CLOSENESS_IN:
			    ret=CLOSENESS_IN_CMD(graph,result);
			    break;
		    case CLOSENESS_OUT:
			    ret=CLOSENESS_OUT_CMD(graph,result);
			    break;
		    case DEGREE_ALL:
			    ret=DEGREE_ALL_CMD(graph,result);
			    break;
		    case DEGREE_IN:
			    ret=DEGREE_IN_CMD(graph,result);
			    break;
		    case DEGREE_OUT:
			    ret=DEGREE_OUT_CMD(graph,result);
			    break;
		    case PAGERANK:
			    ret=PAGERANK_CMD(graph,result);
			    break;
		    default:
			    fprintf(stderr,"%s",USAGE);
			    fprintf(stderr,"Invalid metric\n");
			    exit(EXIT_FAILURE);
	    }	
    }
    if( ret == IGRAPH_ENOMEM){
	    fprintf(stderr,"IGRAPH ERROR\n");
	    exit(EXIT_FAILURE);
    }

    TIMER_STOP;


    ret = igraph_cattribute_VANV(&graph, ID_STRING, igraph_vss_all(),&id_vector); 


    igraph_vector_t extended_results;

    ret = getExtendedResults(&result, &graph, &id_vector, &extended_results);	

    igraph_vector_destroy(&result);

    igraph_vector_destroy(&id_vector);
    unsigned int isolated = 0;
    unsigned int * removed = Calloc(graph.n, sizeof(int));
    unsigned int removed_num = 0;
    unsigned int max_node_to_remove = ((n+100) / 100) * node_perc;

    fout = Fopen(out_file_name, "w+");    
    ftime = timeFileNameGen(out_file_name);
    Free(out_file_name);
    fprintf(ftime,"Search&DeleteBest,Connectivity\n");
 
    fprintf(fout,"%f,%d,%lu\n",(TIMER_ELAPSED /  1000000.0),removed_num,con);
    while ( con > 0 && removed_num < max_node_to_remove){
        
	TIMER_START2; 
        igraph_integer_t id_best = searchIthBestFromId(&graph, &extended_results, removed_num);
        if(id_best == -1){
            fprintf(stderr,"[ERROR]: searchBestByID failed\n");
            exit(1);
        }
        TIMER_STOP;
	removed[removed_num] = VAN(&graph, ID_STRING, id_best);

        ret = igraph_delete_vertices(&graph,igraph_vss_1(id_best));
        if ( ret == IGRAPH_EINVVID){
            fprintf(stderr,"[ERROR]: unable to remove %d\n",VAN(&graph, ID_STRING,id_best));
        }	

        TIMER_STOP2;
        float time_del=(TIMER_ELAPSED2/1000000.0);
	fprintf(fout,"%f,%d,",(TIMER_ELAPSED /  1000000.0), removed[removed_num]);
        removed_num++;
	
 	TIMER_START3;
        con = printNewOutput(removed_num, removed, &isolated, &graph, fout);
        TIMER_STOP3;
        fprintf(ftime,"%f,%f\n",time_del,(TIMER_ELAPSED3/1000000.0));
    }

    igraph_destroy(&graph);
    Free(removed);
    fclose(fout);
    fclose(ftime);
    return 0;
}


