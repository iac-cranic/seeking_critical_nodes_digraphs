//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#ifndef __CRANIC_BF_GRAPH_H_
#define __CRANIC_BF_GRAPH_H_


#include "utility.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>


#ifdef LARGE
	typedef uint64_t vertex_type;
	#define PRIvertex PRIu64
#else
	typedef uint32_t vertex_type;
	#define PRIvertex PRIu32
#endif
	
/* Supported graph file formats */	
#define EDGE_LIST 0
	
/* Supported graph types */
#define DIGRAPH 	0



/*
 * Type: cc_type
 * ------------------------------------------------
 * The CC decomposition of a graph/digraph is implemented using two arrays:
 * 		- <vertex> contains the nodes of the graph;
 *		- <idx> contains the pointers to the nodes in <vertex> (in the same CC).
 *
 * The vertices of CC <i> are in the vector <vertex> from index idx[i] to idx[i+1] excluded. 
 *	
 * Example:
 * 		vertex[idx[i]], ..., vertex[idx[i+1]-1] are in CC <i>.
 *
 *		idx[i+1]-idx[i] is the size of CC <i>.
 *	
 */
typedef struct cc {
    vertex_type *vertex;      /* List of vertices ordered by CC, size |V|=n */
    vertex_type *idx;         /* Indexes of the vertices in the same CC, size <number>+1 */
    vertex_type number;       /* Number of CC */
} cc_type;



/*
 * Type: digraph_type
 * ------------------------------------------------
 * A Digraph is stored using two adjacent lists representing IN and OUT edges.
 * Each list is implemented using 2 arrays:
 * 		- <in_idx> and <in_neighbours> for IN edges;
 *		- <out_idx> and <out_neighbours> for OUT edges.
 * Nodes' ids start from 1.
 *
 * Given a graph G=(V,E) with |V|=n and |E|=m. 
 * The arrays <in_idx> and <out_idx> have size n+2, while the arrays 
 * <in_neighbours> and <out_neighbours> have size m+1. 
 * The positions with index 0 of these 4 arrays are never used for data but only 
 * for implementation purposes.
 * in_idx[0] and out_idx[0] are always set to 0.
 * in_idx[1] and out_idx[1] are always set to 1.
 * The arrays <in_idx> and <out_idx> store the start and end index of the
 * IN and OUT neighbours of each node v of G.
 * Given a node v of G, in_idx[v] is the start index of its IN neighbours
 * in <in_neighbours> and in_idx[v+1] is the end index of its IN neighbours.
 *
 * Example:
 *
 *		out_idx[v+1] - out_idx[v] is the number of OUT neighbours of node v.
 *
 *		out_neighbours[out_idx[v]], out_neighbours[out_idx[v]+1], ..., 
 * 		out_neighbours[out_idx[v+1]-1] are the OUT neighbours of node v.
 *
 * The array <labels> has size n+1 and contains the original ID of the nodes, usually v = labels[v].
 * In case of node removal it could happen that v != labels[v] for some v, because nodes
 * are relabelled after removal.
 */
typedef struct digraph {
	vertex_type	*in_idx;			/* indices of neighbours in <in_neighbours>, size n+2 */
	vertex_type	*in_neighbours;		/* list of IN neighbours (edges), size m+1  */
	vertex_type *out_idx;			/* indices of neighbours in <out_neighbours>, size n+2 */
	vertex_type *out_neighbours;	/* list of OUT neighbours (edges), size m+1  */
    vertex_type *labels;			/* nodes labels, size n+1, labels[i] is the original ID of nodes i */
    vertex_type nodes;				/* number of nodes, |V|=n */
    vertex_type edges;				/* number of edges, |E|=m */
	vertex_type org_nodes_num;		/* number of nodes in the original graph */	
	vertex_type	del_nodes_num;		/* number of deleted nodes */
	vertex_type	*del_nodes_list;    /* list of deleted nodes with their original label */
} digraph_type;



/*
 * Function: get_digraph_from_file
 * ------------------------------------------------
 * The function create a DIGRAPH <digraph_type> from a file. The parameter <file_format>
 * specifies the format of the file. If the file format is not supported a NULL pointer is returned.
 *
 * @param file_name     	The name of the file.
 * @param file_format		The format of the file.
 *
 * @return					The function returns a pointer to a digraph of type digraph_type containing
 * 							the graph. Use the function free_digraph to free the memory. If the file format 
 *                          is not supported or an error occurs a NULL pointer is return      
 */
digraph_type *get_digraph_from_file(const char *file_name, unsigned char file_format);


/*
 * Function: free_digraph
 * ------------------------------------------------
 * The function frees the memory allocated for a digraph. 
 *
 * @param g     	The pointer to the digraph to free.       
 */
void free_digraph(digraph_type *g);


/*
 * Function: get_digraph_scc
 * ------------------------------------------------
 * The function computes the Strongly Connected Components (SCC) of a digraph g. 
 * The function implements the iterative version of the Tarjan algorithm, presented in
 * Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms". 
 *
 * @param g     	The pointer to the digraph. 
 *
 * @return			The funciton returns the pointer to a cc data stracture containing the SCC
 *					composition of the digraph g. 
 */
cc_type *get_digraph_scc(digraph_type *g);

///*
// * Function: print_scc
// * ------------------------------------------------
// * The function prints the CC decomposition of a graph. 
// *
// * @param cc     	 The pointer to the CC to print.  
// */
//void print_cc(cc_type *cc);

/*
 * Function: free_cc
 * ------------------------------------------------
 * The function frees the memory allocated for a cc. 
 *
 * @param cc     	The pointer to the cc to free.       
 */
void free_cc(cc_type *cc);


/*
 * Function: get_pairs_connectivity
 * -----------------------------------
 * Computes the connectivity value of a graph G.
 * Let G a directed or undirected graph, let C_1,C_2,...,C_l be its 
 * (weakly or strongly) connected components.
 * We define the connectivity value of G as:
 *
 *      f(G) = \sum_{i=1}^{l} \binom{C_i}{2}
 *
 * f(G) equals the number of vertex pairs in G that are (weakly or strongly) 
 * connected (i.e., pairwise connectivity value)
 *
 * @param scc       the scc decomposition of the digraph
 *
 * @return          the connectivity value of the digraph
 */
uint64_t get_pairs_connectivity(cc_type *cc);


/*
 * Function: get_digraph_critical_nodes_bruteforce
 * --------------------------------------
 * The function checks which set of nodes of size <k> is critical in the graph <g> and returns
 * a list of them. A set of nodes is critical if it minimaizes the connectivity of <g>.
 * Let G be a digraph, let C_1, C_2, ..., C_l be its strongly connected components.
 * We define the connectivity value of G as:
 *
 *      f(G) = \sum_{i=1}^{l} \binom{C_i}{2}
 *
 * f(G) equals the number of vertex pairs in G that are strongly connected
 * (i.e., pairwise strong connectivity value).
 * For each set of vertecies A of size <k> in V the function creates G\A and computes f(G\A).
 *
 * @param g             	the digraph
 * @param k					the number of verticies to remove 
 * @param max_set			maximum number of critical set to return
 * @param end_conn			residual connectivity after removing any critical set
 * @param tot_sets			total number of critical sets found
 *
 * @return              	An array of critical sets or NULL if <max_set> is set to 0. If an error occurs
 *                          both <end_conn> and <tot_sets> are set to 0 and a NULL pointer is returned.
 */
vertex_type ** get_digraph_critical_nodes_bruteforce ( digraph_type *g, vertex_type k, uint64_t max_set, uint64_t *end_conn, uint64_t *tot_sets ); 


/*
 * Function: delete_nodes_from_digraph
 * --------------------------------------
 * The function removes a set of nodes from a digraph.
 *
 * @param g			the digraph 
 * @param vts		the set of nodes to remove
 * @param n_vts		the number of nodes to remove
 *
 * @return          a copy of the original digraph without the nodes to be removed    	
 */
digraph_type *delete_nodes_from_digraph ( digraph_type *g, vertex_type *vts, vertex_type n_vts);

#endif



