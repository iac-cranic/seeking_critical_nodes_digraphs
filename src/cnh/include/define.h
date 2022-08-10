//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

/* File: define.h */

#ifndef __CNH_DEFINE_H_
#define __CNH_DEFINE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>


#if defined LARGE
typedef unsigned long long VTYPE;
#define  TYPE_FORMAT "%llu"
#else
typedef unsigned int VTYPE;
#define  TYPE_FORMAT "%u"
#endif


#define LTYPE unsigned long long

#define TXT_RED                     "\033[0;31m"
#define TXT_NO_COLOR                "\033[0m"
#define TXT_GREEN                   "\033[0;32m"
#define TXT_BLUE                    "\033[0;34m"
#define TXT_YELLOW                  "\033[1;33m"

#define MAX_BUFFER 4096
#define INIT_POOL_SIZE 64
#define MIN(a,b)  (a < b) ? a : b
#define MAX(a,b)  (a > b) ? a : b

 
/* The vertices of SCC <i> are in vector <vertex> from index <index[i-1]> to <index[i]>,
 * i.e. vertex[index[i-1]], ..., vertex[index[i]] are in SCC <i>.
 * <i> starts from 1, <index[0]> is 0 and size of <index> is <number>+1
 */
typedef struct scc {
    VTYPE *vertex;      /* List of vertices ordered by SCC */
    VTYPE *index;       /* Indexes of the vertices in the same SCC */
    VTYPE number;       /* Number of SCC */
} SCC;

typedef SCC WCC;
typedef SCC CC;


/* A Graph is represented by two matricies in CSR format, one storing IN edges and
 * one storing OUT edges. Each matrix is represented by two arrays:
 * 		- <gin> and <g1in> for IN edges;
 *		- <gout> and <g1out> for OUT edges.
 *
 * Given a graph G=(V,E) with |V|=n and |E|=m, <g1in> and <g1out> have size n+2.
 * They store for each node the indices of nodes with an IN or OUT edge towards them.
 * Their size is n+2 because node ID starts from 1 (0 is not used as node ID).
 * Arrays <gin> and <gout> have size m. They store all nodes having an IN or OUT edge 
 * towards a node in <g1in> and <g1out> respectively.
 * 
 * Example: 
 *	g1out[i+1] - g1out[i] is the number of neighbour of node i (with i!=0).
 *          
 *	gout[g1out[i]], gout[g1out[i]+1], ... , gout[g1out[i+1]-1] are the neighbour
 * 	of node i.
 * 
 *  labels[i] is the original ID of nodes i, usually i = labels[i].
 *  In case of node removal it could happen that i != labels[i] for some i.
 *			
 */
typedef struct Graph {
    VTYPE *gin;				/*  */
    VTYPE *g1in;			/*  */
    VTYPE *gout;			/*  */
    VTYPE *g1out;			/*  */
    VTYPE *labels;			/* Nodes labels, it has size (n+1) */
    VTYPE n;				/* number of nodes, |V|=n */
    VTYPE m;				/* number of edges, |E|=m */
    VTYPE cn; //critical node
    unsigned long long cn_dfo;
    VTYPE *cn_equiv;        /* Array of nodes whose dfo is equal to cn_dfo */
    VTYPE cn_equiv_len;
    VTYPE root;
    VTYPE *removed;
    VTYPE removed_num;
    VTYPE original_n;
    VTYPE analyzed;
    VTYPE randomly_chosen;

} GRAPH;

typedef struct Dfs_Tree {
    // DFS tree arrays
    VTYPE *dfsl;                  // dfs labels in preorder
    VTYPE *ptol;                  // from dfs label to original vertex label
    VTYPE *p;                     // parent array
    VTYPE n;
    VTYPE m;
} DFS_TREE;

typedef struct Clists {
    VTYPE *cycle;                 // list of cycle arcs  (old C_Target)
    VTYPE *last_cycle;            // pointers to the list of cycle arcs (old C_Last)
    VTYPE *cross;                 // cross arcs
    VTYPE *last_cross;
    VTYPE *next;
    VTYPE lastpos;                // index of the last inserted element
    VTYPE n;
    VTYPE m;
} C_EDGE_LISTS;

typedef struct C_Nodes {
    VTYPE *crit_nodes;            // list of critical nodes
    VTYPE len;                    // len of crit_nodes
} C_NODES;

// It must be signed as it may become negative (Is it right??)
typedef struct Bundles {
    LTYPE *sumC;
    LTYPE *sum_sqC;
    LTYPE *sizeC;
} BUNDLES;


// The structure maintains all the subgraphs both already evaluated and to be evaluated
// For each subgraph evaluated (but not selected at the evaluation time), it maintains the subgraph, the critical nodes identified
// along with its dfo_score
typedef struct Graph_Pool {
    VTYPE num_subgraphs;// the number of subgraphs present
    VTYPE start_index;  // starting index of the subgraphs that are not evaluated yet
    VTYPE size;         // the size of subgraphs array
    GRAPH **subgraphs;  // an array of ptr to all subgraphs created by removing critical nodes
    VTYPE idx_best;     // the index of (one) subgraphs (already evaluated) that contains the best critical node according to dfo,
    // w.r.t. to all the other evaluated subgraphs
    VTYPE *idx_best_array;
    VTYPE idx_best_array_len;
    VTYPE isolated_node;
} GRAPH_POOL;

struct GraphStats{
    VTYPE num_isolated_node;
    VTYPE num_subgraphs;
    VTYPE *sizes_subgraphs;
};
typedef struct GraphStats GRAPH_STATS;


/*
 * Function: Malloc
 * --------------------
 * Wrapper of malloc function
 *
 * @param sz: the memory requested size 
 *
 * @return: the pointer to dinamically allocated  memory
 */
void *Malloc(size_t sz);

/*
 * Function: Realloc
 * --------------------
 * Wrapper of realloc function
 *
 * @param ptr: the pointer to memory area to be reallocated
 * @param sz: the memory requested size 
 *
 * @return: the pointer to dinamically allocated  memory
 */
void *Realloc(void *ptr, size_t sz);

/*
 * Function: Free
 * --------------------
 * Wrapper of free function
 *
 * @param ptr: the pointer to the memory to be deallocated 
 *
 * @return
 */

void Free(void *ptr);

/*
 * Function: Calloc
 * --------------------
 * Wrapper of calloc function
 *
 * @param nitems: the number of items of size size to be allocated
 * @param size: the size of each item 
 *
 * @return: the pointer to dinamically allocated  memory
 */
void *Calloc(size_t nitems, size_t size);

/*
 * Function: Strdup
 * --------------------
 * Utility function to manage string copy with dinamically allocated memory
 *
 * @param s: the souce string to be copied
 *  *
 * @return: the pointer to dinamically allocated  memory containing a copy of s
 */
char *Strdup(char *s);

/*
 * Function: Fopen
 * --------------------
 * Wrapper of fopen function
 *
 * @param path: the pathname of the file
 * @param mode: the requested mode 
 *
 * @return: the pointer to the file stream
 */
FILE *Fopen(const char *path, const char *mode);

/*
 * Function: timeFileNameGen
 * --------------------
 * Utility function to open a file name for timing. The timing-filename will be generated
 * starting from the outfilename
 * 
 * @param outfilename: it will be used to generate the timing filename
 *
 *
 * @return: the pointer to file stream
 */
FILE * timeFileNameGen(char *outfilename);

/*
 * Function: usage
 * --------------------
 * Utility function to print the usage message
 * 
 * 
 * @param cmd: the executable name
 *
 *
 * @return
 */
void usage(char *cmd);

#define FATAL_ERROR(message)\
    fprintf(stderr, "[%s][%d] Fatal error: %s\n",__FILE__, __LINE__,  message);\
    exit(111);

#ifdef DEBUG 
#define WARNING(message)\
	fprintf(stderr, "WARNING: %s\n", message);
#else
#define WARNING(message) ;
#endif
	
#define INFO_MESSAGE(message)
/*#define INFO_MESSAGE(message)\
  printTimeStamp();\
    fprintf(stderr, "[INFO]: %s\n", message);\
*/

#endif /* _CNH_DEFINE_H_ */
