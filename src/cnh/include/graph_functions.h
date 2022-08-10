//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

/* FILE: graph_functions.h */

#ifndef __CNH_GRAPH_FUNCTIONS_H_
#define __CNH_GRAPH_FUNCTIONS_H_

#include "define.h"

#ifdef NDEBUG
#define CHECK_GRAPH
#else
#define CHECK_GRAPH(g, func_name)\
    if ((g->m <=0) || (g->n <=0)) {\
        fprintf(stderr, "In function %s", func_name);\
        fprintf(stderr, "Input graph has no elements");\
        exit(111);\
    }
#endif

#ifdef NDEBUG
#define CHECK_STACK_OVERFLOW
#else
#define CHECK_STACK_OVERFLOW\
    if ((2*sid+1) == stacksize) {\
        printf(""TYPE_FORMAT": "TYPE_FORMAT" "TYPE_FORMAT"\n", u, g->g1out[u], g->g1out[u+1]);\
        printf("sid = "TYPE_FORMAT" stacksize = "TYPE_FORMAT"\n", sid, stacksize);\
        fprintf(stderr, "stack overflow!\n");\
        exit(111);\
    }
#endif
	

#define DIMACS 0
#define RMAT 1

/*
 * Function: read_graph
 * --------------------------------------------
 *  It calls the proper function to parse the input file with respect to the fileformat
 *
 *  @param type 	the fileformat, it uses the enum graph_type defined in define.h
 *  @param fname	the input filename
 *  @param n		the number of vertices
 *  @param m		the number of edges
 *  @param root		the root vertex, only for DIMACS format
 *  @param read_edges	if 0 it will only read n and m
 *
 *  @return		the edges list, the number of vertices and edges
 *
 */ 
VTYPE *read_graph(int type, const char *fname, VTYPE *n, VTYPE *m, VTYPE *root, int read_edges);

/*
 * Function: build_graph_datastruct
 * -----------------------------------------
 * It builds the Graph data structure as a CSR.  g1in/out are the row arrays for
 * i, it stores the position in gin/out of the neighbors of
 * each vertices.  gin/out is the column array with neighbors.  n is the number
 * of vertices but labels starts from 1 thus g1in and g1out are allocate with
 * n+2 elements.  the first two elemnts are 0 and the last element is m
 *
 * @param edges		the edges list
 * @param n		the number of vertices
 * @param m 		the number of edges
 * @param root		the root vertex
 *
 * @return 		it returns the CSR representing the input graph
 *
 */
GRAPH *build_graph_datastruct(VTYPE *edges, VTYPE n, VTYPE m, VTYPE root);

/*
 * Function: set_root
 * ------------------------
 * Given the input root vertex r1  and the root vertex r2 read from file, choose one.
 * Always chose r2 if set. If none given
 *
 * @param r1
 * @param r2
 *
 * @return 		the root vertex
 *
 */
VTYPE set_root(VTYPE r1, VTYPE r2);

/*
 * Function: compute_cndp
 * -----------------------------------------
 * It computes the set of k critical nodes accordirg to CNH heuristic proposed in 
 * "Computing Critical Nodes in Directed Graphs. ACM Journal of Experimental Algorihmics 23 (2018)
 * N. Paudel, L. Georgiadis, G. F. Italiano"
 * 
 * @param g		the input graph in CSR format
 * @param k		the number of critical nodes to be removed
 * @param fout	the output file pointer (can be stdout)
 * @param ftime 	the output file pointer (can be stdout)
 *
 * @return		0 on success
 *
 */ 
int compute_cndp(GRAPH *g, VTYPE k, FILE* fout, FILE* ftime);
	
#endif /* GRAPH_FUNCTIONS_H_ */
