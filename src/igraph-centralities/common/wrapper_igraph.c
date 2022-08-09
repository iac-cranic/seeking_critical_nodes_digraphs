//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#include "wrapper_igraph.h"
#include "timer.h"

// read the number of vertices and edges (first line in file)
// and the edge list from the input file;
// stores "nvertices", "nedges" and return the array "edges"
igraph_real_t *read_graph_rmat_igraph(const char *fname, VTYPE *nvertices, VTYPE *nedges, int read_edges) {
	FILE *fp;
	char buffer[BUFFSIZE];
	VTYPE buflen;
	VTYPE comment_count;
	VTYPE nv, ne, i, j;
	VTYPE line_count, n;
	VTYPE nmax;
	igraph_real_t *ed;
	VTYPE self_loops;

	nv = 0;
	ne = 0;
	self_loops = 0;

	fp = Fopen(fname, "r");

	comment_count = 0;
	line_count    = 0;
	n             = 0;
	nmax          = ALLOC_BLOCK; // must be even

	ed = NULL;
	ed = (igraph_real_t *)Malloc(nmax * sizeof(igraph_real_t));


	while (1) {

		// READ LINES
		fgets(buffer, BUFFSIZE, fp);
		buflen = strlen(buffer);
		if (buflen >= BUFFSIZE) {
			fprintf(stderr, "[ERROR]: The line is to long, increase the BUFFSIZE! Exit\n");
			exit(EXIT_FAILURE);
		}

		if (feof(fp)) {
			fclose(fp);
			break;
		} else if (ferror(fp)) {
			fprintf(stderr, "[ERROR]: An error ocurred while reading the file\n");
			fclose(fp);
		}

		// SCAN THE LINE
		if (strstr(buffer, "#") != NULL) {
			if (strstr(buffer, "Nodes:")) {
				sscanf(buffer, "# Nodes: "TYPE_FORMAT" Edges: "TYPE_FORMAT"\n", &nv, &ne);
				if(!read_edges){
					Free(ed);
					return NULL;
				}
			}
			comment_count++;
		} else {
			if ((nv == 0) || (ne == 0)) {
				fprintf(stderr, "[ERROR]: reading the number of vertices or edges in %s\n", fname);
				fprintf(stderr, "[ERROR]: nv = "TYPE_FORMAT" ne = "TYPE_FORMAT"\n", nv, ne);
				exit(EXIT_FAILURE);
			}

			// Read edges
			if (sscanf(buffer, ""TYPE_FORMAT" "TYPE_FORMAT"\n", &i, &j) != 2) {
				fprintf(stderr, "[ERROR]: reading graph arc (%s).\n", fname);
				exit(EXIT_FAILURE);
			}

			if (i >= (nv + 1) || j >= (nv + 1)) {
				fprintf(stderr, "[ERROR]: found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
						fname, (comment_count + line_count), i, j);
				exit(EXIT_FAILURE);
			} else if ((i == 0) || (j == 0)) {
				fprintf(stderr, "[ERROR]: found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
						fname, (comment_count + line_count), i, j);
				fprintf(stderr, "[ERROR]: Label 0 is not allowed, labels must start at 1.");
				exit(EXIT_FAILURE);
			}

			if (n >= nmax) {
				nmax += ALLOC_BLOCK;
				ed   = (igraph_real_t *)Realloc(ed, nmax * sizeof(igraph_real_t));
			}
			if ( i != j) {
				ed[n]   = (double)i-1;
				ed[n + 1] = (double)j-1;
				n       += 2;
			} else {
				self_loops++;
			}
		}
		line_count++;

	}


	if (ne != ((n / 2) + self_loops)) {
		fprintf(stderr, "[ERROR]: reading the input file %s: the number of edges read differ from the number of edges VTYPE the header\n", fname);
		fprintf(stderr, "[ERROR]: nedges header = "TYPE_FORMAT" edges lines read = "TYPE_FORMAT"\n", ne, n / 2);
		exit(EXIT_FAILURE);
	}


	*nedges    = ne;
	*nvertices = nv;

	n = n / 2;
	if (self_loops != 0) {
		fprintf(stderr, "[WARNING]: there were "TYPE_FORMAT" self loops\n", self_loops);
		*nedges = ne - self_loops;
	}
	return ed;
}


// Read flowgraph file type DIMACSGRAPH, source at first line of file followed by the edge list,
// p nv ne src
// a edges
//
// edge list is stored in the array "edges" and function returns the pointer to "edges"
igraph_real_t *read_graph_dimacs_igraph(const char *fname, VTYPE *n, VTYPE *m, VTYPE *root, int read_edges) {
	igraph_real_t  *ledges;                       // The edge list read from file
	VTYPE  u, v, nlines, p, self_loops;
	FILE *fin;
	char line[BUFFSIZE];

	fin = Fopen(fname, "r");

	ledges = NULL;
	nlines = 0;
	p      = 0;
	self_loops = 0;

	while (fgets(line, BUFFSIZE, fin) != NULL) {
		line[strcspn(line, "\r\n")] = 0;

		if (nlines == 0) {
			// Get the number of vertices, arcs and the root vertex: n, m, r
			if (sscanf(line, "p "TYPE_FORMAT" "TYPE_FORMAT" "TYPE_FORMAT"", n, m, root) != 3) {
				fprintf(stderr, "[ERROR]: reading graph size (%s).\n", fname);
				exit(EXIT_FAILURE);
			}
			if(!read_edges){
				return ledges;
			}
			// Allocate memory for the edge list
			ledges = Malloc(2 * (*m) * sizeof(ledges));

			++nlines;
			continue;
		}


		if (sscanf(line, "a "TYPE_FORMAT" "TYPE_FORMAT"", &u, &v) != 2) {
			fprintf(stderr, "[ERROR]: reading graph arc (%s).\n", fname);
			exit(EXIT_FAILURE);
		}

		if (u >= (*n + 1) || v >= (*n + 1)) {
			fprintf(stderr, "[ERROR]: found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
					fname, nlines, u, v);
			exit(EXIT_FAILURE);
		} else if ((u == 0) || (v == 0)) {
			fprintf(stderr, "[ERROR]: found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
					fname, nlines, u, v);
			fprintf(stderr, "[ERROR]: Label 0 is not allowed, labels must start at 1.");
			exit(EXIT_FAILURE);
		}

		if (u == v) {
			self_loops++;
		} else {
			ledges[p++] = (double)u-1;
			ledges[p++] = (double)v-1;
			++nlines;
		}

		if (p > (2 * (*m))) {
			fprintf(stderr, "[ERROR]:  Graph has > "TYPE_FORMAT" arcs.\n", (*m));
			exit(-1);
		}
	}

	if (self_loops > 0) {
		fprintf(stderr, "[WARNING]: there were "TYPE_FORMAT" self loops\n", self_loops);
		*m = p/2;
	}

	fclose(fin);

	return ledges;
}

// uses the enum graph_type defined in define.h
igraph_real_t *read_graph_igraph(int type, const char *fname, VTYPE *n, VTYPE *m, VTYPE *root, int read_edges) {

	igraph_real_t *edges = NULL;

	switch (type) {

		case DIMACS:
			edges =(igraph_real_t*) read_graph_dimacs_igraph(fname, n, m, root, read_edges);
			break;

		case RMAT:
			edges =(igraph_real_t*) read_graph_rmat_igraph(fname, n, m, read_edges);
			break;

		default:
			fprintf(stderr, "[ERROR]:  invalid input graph type %d\n", type);
			exit(EXIT_FAILURE);
			break;

	}

	return edges;
}



#define ELEM_TO_COMPARE_IGRAPH_LIST 2
int list_element_compare_igraph(const void *A, const void *B) {
	int i;
	igraph_real_t *a, *b;

	a = (igraph_real_t*)A;
	b = (igraph_real_t*)B;
	for ( i = 0; i < ELEM_TO_COMPARE_IGRAPH_LIST; i++) {
		if (a[i] < b[i]) {
			//	return -1; //reverting the order (-1 Ascending, 1 Descending)
			return 1;
		} else if ( a[i] > b[i] ) {
			//	return 1; //reverting the order (1 Descending, -1 Ascending) 
			return -1;
		}
	}
	return 0; // It should never reach this point
}

igraph_real_t *sort_list_igraph(igraph_real_t *list, int n_elem) {
	qsort(list, n_elem, ELEM_TO_COMPARE_IGRAPH_LIST * sizeof(igraph_real_t), list_element_compare_igraph);
	return list;
}

igraph_integer_t getBest(igraph_vector_t *result, igraph_t *graph)
{

	VTYPE       i;

	long int result_size = igraph_vector_size(result);
	igraph_real_t *c_result = (igraph_real_t*)Calloc(result_size,sizeof(igraph_real_t));
	igraph_vector_copy_to(result,c_result);
	igraph_real_t *result_handle = (igraph_real_t*)Calloc(2*result_size, sizeof(igraph_real_t));
	for( i = 0 ; i < result_size; i++){
		result_handle[2*i] = c_result[i];
		result_handle[2*i+1] = i;
	} 
	free(c_result);
	sort_list_igraph(result_handle,result_size);
	igraph_integer_t ret = (igraph_integer_t) result_handle[1];
	free(result_handle);
	return ret;
}

unsigned long long getConnectivity(igraph_t *graph){
	igraph_integer_t number_of_cluster = 0;
	igraph_vector_t membership;
	igraph_vector_t csize;
	igraph_vector_init(&membership, graph->n * sizeof(int));
	igraph_vector_init(&csize, graph->n * sizeof(int));
	int ret = igraph_clusters(graph, &membership, &csize, &number_of_cluster, IGRAPH_STRONG);
	if (ret == IGRAPH_EINVAL){
		fprintf(stderr,"[ERROR] internal error\n");
		exit(EXIT_FAILURE);
	}

	int i = 0;
	unsigned long long con = 0;
	igraph_vector_t to_remove;
	igraph_vector_init(&to_remove, graph->n * sizeof(int));
	int to_remove_count = 0;
	for ( i = 0; i < number_of_cluster; i++){
		if( (int)VECTOR(csize)[i] > 1){
			unsigned long long cur = (VECTOR(csize)[i] * (VECTOR(csize)[i] - 1)) / 2;
			con += cur;
		}else {
			int j = 0;
			for(j= 0; j < graph->n; j++){
				if( (int)VECTOR(membership)[j] == i){
					igraph_vector_set(&to_remove, to_remove_count, j);
					to_remove_count++;
				}
			}
		}
	}

	if(to_remove_count > 0 ){    
		if ( ret == IGRAPH_EINVVID){
			fprintf(stderr,"[ERROR]: unable to remove vertices\n");
		}	
	}

	igraph_vector_destroy(&membership);
	igraph_vector_destroy(&csize);
	igraph_vector_destroy(&to_remove);
	return con;
}


unsigned long long printNewOutput(unsigned int num_removed, unsigned int *removed, unsigned int *isolated, igraph_t *graph, FILE *fout){
	unsigned int i = 0;
	unsigned long long connectivity = 0;
	int remove_isolated_vertices = 0;	

	// Computing Connectivity and clusters
	igraph_integer_t number_of_cluster = 0;
	igraph_vector_t membership;
	igraph_vector_t csize;
	igraph_vector_init(&membership, graph->n * sizeof(int));
	igraph_vector_init(&csize, graph->n * sizeof(int));
	int ret = igraph_clusters(graph, &membership, &csize, &number_of_cluster, IGRAPH_STRONG);
	if (ret == IGRAPH_EINVAL){
		fprintf(stderr,"[ERROR] internal error\n");
		exit(EXIT_FAILURE);
	}

	igraph_vector_t to_remove;
	igraph_vector_init(&to_remove, graph->n * sizeof(int));
	int to_remove_count = 0;
        int num_scc=0;
	
	for ( i = 0; i < number_of_cluster; i++){
		if( (int)VECTOR(csize)[i] > 1){
			unsigned long long cur = (VECTOR(csize)[i] * (VECTOR(csize)[i] - 1)) / 2;
			connectivity += cur;
			num_scc++;
		}else {
			int j = 0;
			for(j= 0; j < graph->n; j++){
				if( (int)VECTOR(membership)[j] == i){
					igraph_vector_set(&to_remove, to_remove_count, j);
					to_remove_count++;
				}
			}
		}
	}
	if(to_remove_count > 0){
	    num_scc++;
	}
        fprintf(fout,"%lu,%u,%u,%d",connectivity,graph->n,igraph_ecount(graph),num_scc);	
	//for ( i = 0; i < number_of_cluster; i++){
	//	if( (int)VECTOR(csize)[i] > 1){
	//		fprintf(fout,",%d",(int)VECTOR(csize)[i]);
	//	}
	//}
	fprintf(fout,",%d\n",to_remove_count);

	*isolated += to_remove_count;

	if(to_remove_count > 0  && remove_isolated_vertices){    
		ret = igraph_delete_vertices(graph,igraph_vss_vector(&to_remove));
		if ( ret == IGRAPH_EINVVID){
			fprintf(stderr,"[ERROR]: unable to remove vertices\n");
			exit(EXIT_FAILURE);
		}	
	}
	igraph_vector_destroy(&membership);
	igraph_vector_destroy(&csize);
	igraph_vector_destroy(&to_remove);

	return connectivity;
}


int getExtendedResults(igraph_vector_t *res, igraph_t *graph, igraph_vector_t *id, igraph_vector_t *ex_res){


	long int result_size = igraph_vector_size(res);
	int i = 0;

	igraph_real_t *c_result = (igraph_real_t*)Calloc(result_size,sizeof(igraph_real_t));
	igraph_vector_copy_to(res,c_result);
	igraph_real_t *result_handle = (igraph_real_t*)Calloc(2*result_size, sizeof(igraph_real_t));
	for( i = 0 ; i < result_size; i++){
		result_handle[2*i] = c_result[i];
		result_handle[2*i+1] = (int)VAN(graph, ID_STRING, i);
	} 
	free(c_result);
	sort_list_igraph(result_handle,result_size);

	int ret = igraph_vector_init_copy(ex_res, result_handle, 2*result_size);
	if(ret == IGRAPH_ENOMEM){
		fprintf(stderr,"[ERROR]: igraph_vector_init_copy\n");
		return -1;
	}
	free(result_handle);
	return 0;

}

igraph_integer_t searchIthBestFromId(igraph_t *graph, igraph_vector_t *ex_res, int idx){


	int i = 0;
	int x = (int)VECTOR(*ex_res)[2*idx + 1];
	for(i = 0; i < graph->n; i++){
		if(x == (int)VAN(graph, ID_STRING,i)){
			return i;
		}
	}
	return -1;
}
