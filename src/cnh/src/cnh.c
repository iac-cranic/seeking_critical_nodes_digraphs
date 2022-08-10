//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

/* File: cnh.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include "define.h"
#include "timer.h"
#include "graph_functions.h"



// GLOBAL VARS
// Dominator arrays
int *d;       // dominators

// Loop Nesting Tree
int *h;       // headers

// Type of graph files supported
const char *GTYPE_NAMES[2] = {"DIMACS", "RMAT"};


unsigned  long int f(VTYPE n) {
    unsigned long int a = (unsigned  long int)(n) * (n - 1);
    a = a / 2;
    return a;
}

// Goals
// 1) Report all the strongly connected components obtained after deleting any single vertex
// 2) Compute the number of strongly connected components obtained after deleted a single vertex for all vertices
// 3) Compute the size of the largest and of the smallest strongly connected component obtained after deleting a single vertex
// 	  for all vertex
int main(int argc, char **argv) {

    int         gtype;                  // Graph type (see enum graph_type)
    char        c;                      // getopt char
    int		read_edges;				// Determine if read_graph should read also the edges list

    VTYPE        *edges;                // edge list read from file
    GRAPH        *g;                    // Input Graph

    VTYPE        root;                  // Source vertex of the flow graph
    VTYPE        read_root;             // Soource vertex read from file (DIMACS)
    VTYPE        n, m;                  // n -> number of vertices, m -> number of edges
    int 	 k;      		// number of vertices to remove
    

    // Others
    char *graph_file_name;        // input file name
    char *out_file_name;	  // output file name

    FILE *fout;
    FILE *ftime;

    // Initialize
    gtype     = -1;
    read_edges=  1;
    n         =  0;
    m         =  0;
    root      = -1;
    read_root = -1;
    k         =  0;

    g        = NULL;

    graph_file_name = NULL;
    out_file_name = NULL;
    
    fout = stdout;
    ftime = stdout;

    // Parse command line
    /* Manage command line */
    if (argc < 4 || argc > 9) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    while ((c = getopt(argc, argv, "f:t:r:k:o:h")) != EOF) {
        switch (c) {
            // Graph file name
            case 'f':
                graph_file_name = Strdup(optarg);
                break;
                // Graph file type
            case 't':
                gtype = atoi(optarg);
                break;
                // Root vertex
            case 'r':
                root = atoi(optarg);
                break;
                // Output file name 
            case 'o':
                out_file_name = Strdup(optarg);
                break;
                // Number of vertex to remove
            case 'k':
                k = atoi(optarg);
                break;
                // Help
            case 'h':       // Help
                usage(argv[0]);
                exit(EXIT_FAILURE);
                break;
            default:
                usage(argv[0]);
                exit(EXIT_FAILURE);
                break;
        }
    }

    if (gtype) {
        if (graph_file_name == NULL) {
            fprintf(stderr, "option -t require -f <input file>\n");
            exit(EXIT_FAILURE);
        }
    }
    if (graph_file_name) {
        if (gtype < 0) {
            fprintf(stderr, "option -f require -t <graph type>\n");
            exit(EXIT_FAILURE);
        }
    }
    if (gtype > RMAT) {
        fprintf(stderr, "file type not recognized!\n");
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    if (k < 0) {
        fprintf(stderr, "k must be greaten or equal  than 0\n");
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    if (!out_file_name){
	fprintf(stderr,"the output will be printed on stdout\n");
    } else {
	fout = Fopen(out_file_name,"w+");
	ftime=timeFileNameGen(out_file_name);
    }
    printf("input file name: %s\n", graph_file_name);
    printf("graph type: %d: %s\n", gtype, GTYPE_NAMES[gtype]);


    // Read the graph from file and store the edges in an edge list
    edges = read_graph(gtype, graph_file_name, &n, &m, &read_root, read_edges);
    root  = set_root(root, read_root);
    
    // Build a CSR data structure to represent the graph
    g = build_graph_datastruct(edges, n, m, root);
    free(edges);
    compute_cndp(g,k,fout,ftime);

    free(out_file_name);
    free(graph_file_name);
    return 0;
}

