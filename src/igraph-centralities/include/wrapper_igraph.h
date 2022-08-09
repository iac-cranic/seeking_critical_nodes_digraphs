//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#ifndef __IGRAPH_CENTRALITIES_WRAPPER_IGRAPH_H_
#define __IGRAPH_CENTRALITIES_WRAPPER_IGRAPH_H_

#include <igraph.h>
#include "define.h"

#define ID_STRING "id"

igraph_real_t *read_graph_rmat_igraph(const char *fname, VTYPE *nvertices, VTYPE *nedges, int read_edges);
igraph_real_t *read_graph_dimacs_igraph(const char *fname, VTYPE *n, VTYPE *m, VTYPE *root, int read_edges);
igraph_real_t *read_graph_igraph(int type, const char *fname, VTYPE *n, VTYPE *m, VTYPE *root, int read_edges);

int list_element_compare_igraph(const void *A, const void *B);
igraph_real_t *sort_list_igraph(igraph_real_t *list, int n_elem);
igraph_integer_t getBest(igraph_vector_t *result, igraph_t *graph);
unsigned long long getConnectivity(igraph_t *graph);
unsigned long long printNewOutput(unsigned int num_removed, unsigned int *removed, unsigned int *isolated, igraph_t *graph, FILE *fout);
int getExtendedResults(igraph_vector_t *res, igraph_t *graph, igraph_vector_t *id, igraph_vector_t *ex_res);
igraph_integer_t searchIthBestFromId(igraph_t *graph, igraph_vector_t *ex_res, int idx);
#define ALLOC_BLOCK (2*1024)
#define BUFFSIZE     1024

#endif
