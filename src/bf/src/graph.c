//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#include "graph.h"


#define BUFFSIZE     1024
/*
 * Function: get_edge_list_from_file
 * ------------------------------------------------
 * The function read a graph file in an edge list format. It returns an array containing
 * the list of edges of the graph read from the file. The function also fill the two parameters
 * <n_nodes> and <n_edges> with the number of nodes and edges of the graph respectivally.
 *
 * @param file_name		The name of the file
 * @param n_nodes     	The parameter will contains the number of nodes in the graph
 * @param n_edges     	The parameter will contains the number of edges in the graph
 *
 * @return          	The function allocates and returns a pointer to an array containing the 
 *						list of edges of the graph. Returns a NULL pointer if an error occurs.
 */
static vertex_type *get_edge_list_from_file(const char *file_name, vertex_type *n_nodes, vertex_type *n_edges){
	FILE *fp=NULL;
	char buffer[BUFFSIZE];
	vertex_type *edge_list=NULL;
	vertex_type edge_idx=0;
	vertex_type self_loops=0;
	int err=0, n_line=0, n_comments=0;
	vertex_type nodes=0; 
	vertex_type edges=0;
	
	fp = Fopen(file_name, "r", __FILE__, __LINE__, "Opening an edge list file.");
	
	while ( fgets(buffer, BUFFSIZE, fp)){ /* READ a LINE */
		int buflen = 0;
		buflen = strlen(buffer);
        if ( buffer[buflen-1] != '\n' && !feof(fp) ) { // Check that the buffer is big enough to read a line
            fprintf(stderr, "(%s)[ERROR] File %s. The line is too long, increase the BUFFSIZE! Exit\n",get_time_string(), file_name);
			fclose(fp);
			return NULL;
        }
		n_line++;
		
		if (strchr(buffer, '#') != NULL) { // The line it's a comment
			n_comments++;
            if (strstr(buffer, "Nodes:") ) { // The line contains the number of Nodes and Edges
                err = sscanf(buffer, "# Nodes: %"PRIvertex" Edges: %"PRIvertex"\n", &nodes, &edges);
				if (err != 2 || err == EOF){
					fprintf(stderr, "(%s)[ERROR] Error reading file %s.\n",get_time_string(), file_name);
					fclose(fp);
					return NULL;
				} /* Allocate the memory for the list of edges */
				edge_list = (vertex_type*) Malloc(sizeof(vertex_type)*edges*2,__FILE__, __LINE__, "Allocating the edge list.");
            }
		}else{ // Read a new edge
			vertex_type i, j;
			if( nodes == 0 || edges == 0){
				fprintf(stderr, "(%s)[ERROR] File %s. You must specify the number of Nodes and Edges in the graph before listing the edges.\n",get_time_string(), file_name);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			sscanf(buffer, "%"PRIvertex" %"PRIvertex"\n", &i, &j);
			if (err != 2 || err == EOF){
				fprintf(stderr, "(%s)[ERROR] Error reading file %s.\n", get_time_string(), file_name);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			if (i > nodes || j > nodes){
				fprintf(stderr, "(%s)[ERROR] File %s. Vertex id bigger than the number of vertices in the graph: edge %"PRIvertex" %"PRIvertex", number of nodes %"PRIvertex".\n",get_time_string(), file_name, i, j, nodes);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			if (i == 0 || j == 0){
				fprintf(stderr, "(%s)[ERROR] File %s. Vertex id 0 not allowed, vertex id must start from 1.\n",get_time_string(), file_name);
				fclose(fp);
				free(edge_list);
				return NULL;
			}
			
			if ( i != j) { // Ignore self-loops
				if( edge_idx == (edges*2) ){
					fprintf(stderr, "(%s)[ERROR] File %s. Too many edges listed in the file (%"PRIvertex"), only %"PRIvertex" were expected\n",get_time_string(), file_name, edge_idx/2, edges );
					fclose(fp);
					free(edge_list);
					return NULL;
				}
				edge_list[edge_idx]=i;
				edge_list[edge_idx+1]=j;
				edge_idx +=2;
			} else { self_loops++; }
		}
	}	
 
    if (ferror(fp)){ /* Check reading errors */
		perror("Reading error");
		fclose(fp);
		return NULL;
	}
		
	if( ((edge_idx/2)+self_loops) > edges ){
		fprintf(stderr, "(%s)[ERROR] File %s. Too many edges listed in the file, only %"PRIvertex" were expected\n",get_time_string(), file_name, edges );
		fclose(fp);
		free(edge_list);
		return NULL;
	}
	fclose(fp);
	
	if ( ((edge_idx/2)+self_loops) < edges ){
		fprintf(stderr, "(%s)[ERROR] File %s. Not enough edges listed in the file, %"PRIvertex" were expected\n",get_time_string(), file_name, edges );
		free(edge_list);
		return NULL;
	} 
	
	*n_edges=edges-self_loops;
	*n_nodes=nodes;
	
	return edge_list;
}




/*
 * Function: get_digraph_from_edge_list
 * ------------------------------------------------
 * The function create a DIGRAPH <digraph_type> from an edge list. The edge list is an array
 * containing the edges of the graph, position i and i+1 are the src and dst of each edge.
 *
 * @param edge_list     	The array containing the list of edges.
 * @param n_nodes			Number of nodes in the graph.
 * @param n_edges			Number of edges in the graph.
 *
 * @return					The function returns a pointer to a digraph of type digraph_type containing
 * 							the graph. Use the function free_digraph to free the memory.       
 */
static digraph_type *get_digraph_from_edge_list(vertex_type *edge_list, vertex_type n_nodes, vertex_type n_edges){
	vertex_type	*in_idx;			
	vertex_type	*in_neighbours;		
	vertex_type *out_idx;			
	vertex_type *out_neighbours;	
	vertex_type *labels;
	digraph_type *g;
	vertex_type *next_in;
	vertex_type *next_out;
	
	// Allocate the memory for the graph
	g = (digraph_type *) Malloc(sizeof(digraph_type),__FILE__, __LINE__, "Graph allocation.");
	g->nodes = n_nodes;
	g->edges = n_edges;
	labels = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "Graph labels.");
	in_idx = (vertex_type *) Calloc((n_nodes+2),sizeof(vertex_type),__FILE__, __LINE__, "Graph in_idx.");
	out_idx = (vertex_type *) Calloc((n_nodes+2),sizeof(vertex_type),__FILE__, __LINE__, "Graph out_idx.");
	in_neighbours = (vertex_type *) Malloc(sizeof(vertex_type)*(n_edges+1),__FILE__, __LINE__, "Graph in_neighbours.");
	out_neighbours = (vertex_type *) Malloc(sizeof(vertex_type)*(n_edges+1),__FILE__, __LINE__, "Graph out_neighbours.");
	// Allocate the memory for temporary structures
	next_in = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "next_in");
	next_out = (vertex_type *) Calloc(n_nodes+1,sizeof(vertex_type),__FILE__, __LINE__, "next_out");
	
	in_idx[1] = out_idx[1] = 1;
	in_neighbours[0] = out_neighbours[0] = 0;
	
	// We write only indexes between 2 and n+1 of arrays <in_idx> and <out_idx>.
	for (vertex_type i = 0; i < n_edges*2; i=i+2) { // Create in_idx and out_idx arrays 
		in_idx[ edge_list[i+1] + 1 ]++;
		out_idx[ edge_list[i] + 1 ]++;
	}
	for (vertex_type i = 0; i < n_nodes+1; i++) { // Create in_idx and out_idx arrays 
		in_idx[i+1] += in_idx[i];
		out_idx[i+1] += out_idx[i];
		next_in[i] = in_idx[i];
		next_out[i] = out_idx[i];
		labels[i] = i;
	}

	for (vertex_type i = 0; i < n_edges*2; i=i+2) { // Create in_neighbours and out_neighbours arrays
		in_neighbours[ next_in[edge_list[i+1]]++ ] = edge_list[i];
		out_neighbours[ next_out[ edge_list[i]]++ ] = edge_list[i+1];
	}
	free(next_in);
	free(next_out);
	
	g->labels = labels;
	g->in_idx = in_idx;
	g->out_idx = out_idx;
	g->in_neighbours = in_neighbours;
	g->out_neighbours = out_neighbours;
    g->del_nodes_list = NULL;
    g->del_nodes_num = 0;
    g->org_nodes_num = g->nodes;
	
	return g;
}



/*
 * Function: compare_edge
 * ------------------------------------------------     
 * Compare two edges based on the ID of their verticies. 
 *
 * @param edge1   First edge to compare
 * @param edge2   Second edge to compare
 *
 * @return		0 if the edges are equal, 1 if edge1 is greater than edge2 
 *              and -1 if edge2 is greater than edge1
 */
static int compare_edge(const void *edge1, const void *edge2) {
    vertex_type *ed1 = (vertex_type *) edge1;
    vertex_type *ed2 = (vertex_type *) edge2;

    if (ed1[0] < ed2[0]) return -1;
    if (ed1[0] > ed2[0]) return  1;
    if (ed1[1] < ed2[1]) return -1;
    if (ed1[1] > ed2[1]) return  1;
    return 0;
}
/*
 * Function: remove_multiple_edges
 * ------------------------------------------------     
 * The function removes multiple edges, keeping only one in case of duplicates
 * between nodes.
 *
 * @param edge_list   An array containing the list of edges, each edge is stored 
 *                    in two adjacent position of the array
 * @param n_edges     Number of edges
 *
 * @return		Number of edges kept after removing duplicates.
 */
static vertex_type remove_multiple_edges(vertex_type *edge_list, vertex_type n_edges){
	vertex_type kept, i;
		
	kept=0;	
	qsort(edge_list, n_edges, sizeof(vertex_type[2]), compare_edge);
	
    for( i=0; i < (n_edges-1)*2; i=i+2){
		if ( edge_list[i]==edge_list[i+2] && edge_list[i+1]==edge_list[i+3]){ continue; }		
        //fprintf(stderr,"%"PRIvertex" - %"PRIvertex"\n", i, kept);
		edge_list[kept] = edge_list[i];
		edge_list[kept+1] = edge_list[i+1];
		kept = kept+2;
	}
    //fprintf(stderr,"%"PRIvertex" - %"PRIvertex"\n", i, kept);
	edge_list[kept] = edge_list[i];
	edge_list[kept+1] = edge_list[i+1];
	kept = kept+2;
	
	fprintf(stderr, "(%s)[INFO] Number of multiple edges removed %"PRIvertex"\n", get_time_string(), n_edges-(kept/2));
	return kept/2;
}

/*
 * Function: get_digraph_from_file
 * ------------------------------------------------
 * The function create a DIGRAPH <digraph_type> from a file. The parameter <file_format>
 * specifies the format of the file. If the file format is not supported a NULL pointer is returned
 *
 * @param file_name     	The name of the file.
 * @param file_format		The format of the file.
 *
 * @return					The function returns a pointer to a digraph of type digraph_type containing
 * 							the graph. Use the function free_digraph to free the memory. If the file format 
 *                          is not supported or an error occurs a NULL pointer is returned      
 */
digraph_type *get_digraph_from_file(const char *file_name, unsigned char file_format){
	vertex_type n_nodes; 
	vertex_type n_edges;
	vertex_type *edge_list;
	digraph_type *g;
	
	fprintf(stderr, "(%s)[INFO] Opening file %s\n",get_time_string(), file_name);
	switch (file_format){
		case EDGE_LIST:
			edge_list = get_edge_list_from_file(file_name, &n_nodes, &n_edges);
			break;
		default:
			fprintf(stderr, "(%s)[ERROR] Unknown input file format %u\n",get_time_string(), file_format);
			return NULL;
	}
	if (edge_list == NULL){ /* An error occurred while reading the file */
		return NULL;
	}
	fprintf(stderr, "(%s)[INFO] File %s. The graph has %"PRIvertex" Nodes and %"PRIvertex" Edges\n",get_time_string(), file_name, n_nodes, n_edges);
	n_edges = remove_multiple_edges(edge_list, n_edges);
	fprintf(stderr, "(%s)[INFO] File %s. The graph has %"PRIvertex" Nodes and %"PRIvertex" Edges\n",get_time_string(), file_name, n_nodes, n_edges);
	g = get_digraph_from_edge_list(edge_list, n_nodes, n_edges);
	free(edge_list);
	fprintf(stderr, "(%s)[DEBUG] Graph building completed. \n",get_time_string());
	
	return g;
}


/*
 * Function: free_digraph
 * ------------------------------------------------
 * The function frees the memory allocated for a digraph. 
 *
 * @param g     	The pointer to the digraph to free.       
 */
void free_digraph(digraph_type *g){
	if (g != NULL){
		free(g->labels);
		free(g->in_idx);
		free(g->out_idx);
		free(g->in_neighbours);
		free(g->out_neighbours);
		free(g->del_nodes_list);
		free(g);
	}
	return;
}

/*
 * Function: get_digraph_scc
 * ------------------------------------------------
 * The function computes the Strongly Connected Components (SCC) of a digraph g. 
 * The function implements the iterative version of the Tarjan algorithm, presented in
 * Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms". 
 *
 * @param g     	The pointer to the digraph. 
 *
 * @return			The function returns the pointer to a cc data stracture containing the SCC
 *					composition of the digraph g. 
 */
cc_type *get_digraph_scc(digraph_type *g){
	vertex_type *lowlink, *number;	
	vertex_type *scc_stack, *stack;	
	char *node_in_scc_stack;
	cc_type *scc;	
	vertex_type stack_idx, last_unvisited, dfs_idx, scc_stack_idx, n_nodes;
	char proceed;
	
	n_nodes = g->nodes; 
	lowlink = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "Lowlink array.");
	number = (vertex_type *) Calloc( n_nodes+1, sizeof(vertex_type), __FILE__, __LINE__, "Number array.");
	stack = (vertex_type *) Malloc(sizeof(vertex_type)*(g->edges*2),__FILE__, __LINE__, "Stack array.");
	scc_stack = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1),__FILE__, __LINE__, "SCC Stack array.");
	node_in_scc_stack = (char *) Calloc((n_nodes+1), sizeof(char), __FILE__, __LINE__, "Check SCC Stack array.");
	scc = (cc_type*) Malloc(sizeof(cc_type), __FILE__, __LINE__, "SCC." );
	scc->vertex = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "scc->vertex array.");
	scc->idx = (vertex_type *) Malloc(sizeof(vertex_type)*(n_nodes+1), __FILE__, __LINE__, "scc->idx array.");
	
	proceed = 1;
	last_unvisited = 1;
	stack_idx = 1;
	dfs_idx = 1;
	scc_stack_idx = 0;
	stack[stack_idx] = 1;
	scc->number = 0;
	scc->idx[0] = 0;
	
	while(proceed == 1){
		vertex_type v;
		vertex_type w;
		v = stack[stack_idx];
		if(number[v] == 0){ /* _ New node to visit, we keep it in the stack _ */
			scc_stack_idx++;	/* put the node in the stack of points */
			scc_stack[scc_stack_idx]=v;
			node_in_scc_stack[v]=1;
			lowlink[v]=dfs_idx;	/* set lowlink and number of the node */
			number[v]=dfs_idx;
			dfs_idx++;
			for(vertex_type j = g->out_idx[v+1]-1; j >= g->out_idx[v]; j--){ /* For each neighbour of the node */
				w = g->out_neighbours[j];
				if( number[w] == 0 ){ /* Put in the stack the nodes we didn't visit yet */
					stack_idx++;
					stack[stack_idx] = w;	
				}else if ( number[w] < number[v] && node_in_scc_stack[w] == 1){  /* Update the lowlink of v if the neighbour has been visited and it is in the stack of points  */
					lowlink[v] = MIN( lowlink[v], number[w]);
				}
			}
		}else if(node_in_scc_stack[v] == 1) { /* _ A node we visited and it is still in the stack of points _ */
			stack_idx--;	/* Remove the node from the stack */
			for(vertex_type j = g->out_idx[v+1]-1; j >= g->out_idx[v]; j--){ /* For each neighbour of the node */
				w = g->out_neighbours[j];
				if ( node_in_scc_stack[w] == 1 ){ /* Update the lowlink of v if the neighbour w has been visited and w is in the stack of points  */
					lowlink[v] = MIN( lowlink[v], lowlink[w]);
				}
			}
			if(lowlink[v] == number[v] ){ /* v is the root of a tree, v is the root of a SCC */
				vertex_type scc_size = 0;
				vertex_type scc_idx = scc->idx[scc->number];
				scc->number++;
				w = scc_stack[scc_stack_idx];
				while( number[w] >= number[v] ){
					scc_size++;
					node_in_scc_stack[w] = 0;
					scc->vertex[scc_idx] = w; 	/* update scc data structure */
					scc_idx++;
					if ( scc_stack_idx == 1 ) { break; }
					scc_stack_idx--; 	/* read next vertex that is possibly in the scc */
					w = scc_stack[scc_stack_idx];
				}	
				scc->idx[scc->number] = scc_idx;
			}
		}else{ /* ___ A node we visited and it is NO MORE in the stack of points (scc_stack) ___ */
			stack_idx--;	/* Remove the node from the stack */
		}		
		if( stack_idx == 0 ){ /* This was the last node in the stack */
			if ( dfs_idx == n_nodes+1 ){ /* We visited all nodes */
				proceed = 0; 
			}else{
				for(vertex_type i = last_unvisited; i < n_nodes+1; i++ ){  /* There are other nodes to visit */
					if(number[i] == 0){
						stack_idx++;
						stack[stack_idx] = i;
						last_unvisited = i+1;
						break;
					}
				}
			}
		}	
	}
	scc->idx = (vertex_type*) Realloc(scc->idx, sizeof(vertex_type)*(scc->number+1), __FILE__, __LINE__, "scc->idx array.");
    free(stack);
    free(scc_stack);
    free(lowlink);
    free(number);
    free(node_in_scc_stack);
	return scc;	
}

/*
 * Function: free_cc
 * ------------------------------------------------
 * The function frees the memory allocated for a CC. 
 *
 * @param cc     	The pointer to the CC to free.       
 */
void free_cc (cc_type *cc){
	if(cc != NULL){
		free(cc->vertex);
		free(cc->idx);
		free(cc);
	}
}


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
uint64_t get_pairs_connectivity(cc_type *cc) {
    uint64_t conn_val;		  /* the connectivity value */

    conn_val = 0;
    for ( vertex_type i = 1; i <= cc->number; i++) {
		uint64_t c_size = 0;
		c_size = cc->idx[i] - cc->idx[i - 1];
        conn_val += c_size*(c_size-1)/2 ;
    }
    return conn_val;
}







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
digraph_type *delete_nodes_from_digraph(digraph_type *g, vertex_type *vts, vertex_type n_vts) {
    digraph_type *g_minus;				/* The new graph G\A */
	vertex_type n_nodes, n_edges, n_kept_nodes;
	char *to_remove;
	vertex_type *new_id;

	n_nodes = g->nodes;
	n_edges = g->edges;
	n_kept_nodes = n_nodes - n_vts;
    g_minus = (digraph_type *) Calloc( 1, sizeof(digraph_type),__FILE__, __LINE__, "Calloc of g_minus");
    g_minus->nodes = n_kept_nodes;
	g_minus->edges = n_edges;
    g_minus->labels = (vertex_type *) Malloc( (n_kept_nodes + 1)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->labels");
    g_minus->in_idx = (vertex_type *) Malloc( (n_kept_nodes + 2)*sizeof(vertex_type), __FILE__, __LINE__,"Malloc of g_minus->in_idx");
    g_minus->out_idx = (vertex_type *) Malloc( (n_kept_nodes + 2)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->out_idx");
    g_minus->in_neighbours = (vertex_type *) Malloc( (n_edges + 1)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->in_neighbours");
    g_minus->out_neighbours = (vertex_type *) Malloc( (n_edges + 1)*sizeof(vertex_type), __FILE__, __LINE__,"Malloc of g_minus->out_neighbours");
    g_minus->del_nodes_num = g->del_nodes_num + n_vts;
    g_minus->del_nodes_list = (vertex_type *) Malloc( (g_minus->del_nodes_num)*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of g_minus->del_nodes_list");
    g_minus->org_nodes_num = g->org_nodes_num;
	to_remove = (char *) Calloc(n_nodes+1, sizeof(char), __FILE__, __LINE__,"Calloc to_remove array");
	
	if (g->del_nodes_num > 0 ){
		memcpy(g_minus->del_nodes_list, g->del_nodes_list, (g->del_nodes_num)*sizeof(vertex_type));
    	for(vertex_type x = g->del_nodes_num, y=0; x < g_minus->del_nodes_num; x++, y++){
			g_minus->del_nodes_list[x] = g->labels[vts[y]];
    	}
	}else{
		memcpy(g_minus->del_nodes_list, vts, n_vts*sizeof(vertex_type));
	}
	
    for(vertex_type x = 0; x < n_vts; x++){
		if(vts[x] > n_nodes){
			return NULL;
		}
		to_remove[vts[x]] = 1;
    }
	new_id = (vertex_type *) Calloc(n_nodes+1, sizeof(vertex_type),__FILE__, __LINE__, "Calloc to_remove array");
	
	/* Remove out edges */
	vertex_type j, h, del_edges, kept_edges;
    for (j = 1, kept_edges = 1, h = 1, del_edges = 0; j < n_nodes+1; j++) { 
		/* Checks all vertices, scans out_idx and out_neighbours */
        if ( to_remove[j] == 1) { /* We remove <j> and all edges starting from <j> */
			del_edges += g->out_idx[j+1] - g->out_idx[j];
        }else{   /* We do not remove <j>, <h> is the ID of <j> in the new graph. */
            g_minus->labels[h] = g->labels[j];
            g_minus->out_idx[h] = g->out_idx[j] - del_edges;
			new_id[j] = h;
            h++;
			/* We do not remove <j>, we checks if we can copy its neighbors in the new graph */
            for (vertex_type i = g->out_idx[j]; i < g->out_idx[j + 1]; i++) { /* Scan all neighbors <u> of <j> */
				vertex_type u;
                u = g->out_neighbours[i];
				if (to_remove[u] == 0){ /* we keep the neighbour */
					g_minus->out_neighbours[kept_edges] = u;
					kept_edges++;
				}else{ del_edges++; }
            }
        }
    }
    g_minus->out_idx[h] = g->out_idx[j] - del_edges; /* Update out_idx last index [n+1]*/
    for(vertex_type x = 1; x < kept_edges; x++){ /* update nodes id */
		g_minus->out_neighbours[x] = new_id[g_minus->out_neighbours[x]];
    }
	
	/* Remove in edges */
    for (j = 1, kept_edges = 1, h = 1, del_edges = 0; j < n_nodes+1; j++) { /* Checks all vertices */
        if ( to_remove[j] == 1) { /* We remove <j> and all edges starting from <j> */
            del_edges += g->in_idx[j + 1] - g->in_idx[j];
        } else {  /* We do not remove <j>, <h> is the ID of <j> in the new graph. */
            g_minus->in_idx[h] = g->in_idx[j] - del_edges;
            h++;
            for (vertex_type i = g->in_idx[j]; i < g->in_idx[j + 1]; i++) { /* Scan all neighbors <u> of <j> */
				vertex_type u;
                u = g->in_neighbours[i];
				if (to_remove[u] == 0){ /* we keep the neighbour */
					g_minus->in_neighbours[kept_edges] = new_id[u];
					kept_edges++;
				}else{ del_edges++; }
            }
        }
    }
    g_minus->in_idx[h] = g->in_idx[j] - del_edges; /* Update in_idx last index [n+1]*/

    g_minus->edges -= del_edges;
    g_minus->in_neighbours = Realloc( g_minus->in_neighbours, (kept_edges) * sizeof(vertex_type), __FILE__, __LINE__, "Realloc of g_minus->in_neighbours" );
    g_minus->out_neighbours = Realloc( g_minus->out_neighbours,  (kept_edges) * sizeof(vertex_type), __FILE__, __LINE__, "Realloc of g_minus->out_neighbours" );
	free(to_remove);
    free(new_id);
	
    return g_minus;
}



/*
 * Function: get_digraph_critical_nodes_bruteforce
 * --------------------------------------
 * The function checks which set of nodes of size <k> is critical in the graph <g> and returns
 * a list of them. A set of nodes is critical if it minimaizes the connectivity of <g>.
 * Let G a digraph, let C_1, C_2, ..., C_l be its strongly connected components.
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
vertex_type ** get_digraph_critical_nodes_bruteforce ( digraph_type *g, vertex_type k, uint64_t max_set, uint64_t *end_conn, uint64_t *tot_sets ) {
	vertex_type n;					/* Number of verticies in G 								*/
	vertex_type **opt_set;          /* List of optimal critical sets found 						*/
    uint64_t n_crtl_set;            /* Number of critical sets found 							*/
    uint64_t crtl_set_idx;          /* Index of the next critical set to insert in <opt_set>	*/    

    n_crtl_set = 0;
    crtl_set_idx = 0;
	n = g->nodes;
	if( k > n ){ *tot_sets = 0; *end_conn=0; return NULL; }
	if (max_set != 0){ // number of optimal set to store
    	opt_set = (vertex_type **) Malloc( max_set*sizeof(vertex_type *),__FILE__, __LINE__, "Malloc of opt_set" );
    	for (vertex_type i = 0; i < max_set; i++){
        	opt_set[i] = (vertex_type *) Malloc( k*sizeof(vertex_type),__FILE__, __LINE__, "Malloc of opt_set[i]" );
    	}
	}else{ opt_set = NULL; }

    /* Compute initial values SCC and Connectivity */
    cc_type *scc;                          /* Number of scc found in G */
	uint64_t min_conn;                      /* Min connectivity value for set of size <k> */
    vertex_type *set_gen, l, j, lst_e;
    digraph_type *diff_g;                   /* G after verticies removal */
    uint64_t conn, idx_s;
	
    /* Start Generating the subsets */
	set_gen = Malloc( k*sizeof(vertex_type), __FILE__, __LINE__,"Malloc of set_gen" );
	for ( j = 0; j < k; j++){
		set_gen[j] = j+1;
	}
	l = k-1;
	lst_e = k-1;
	idx_s = 0;		// Try to remove the first subset
    diff_g = delete_nodes_from_digraph( g, set_gen, k );
    scc = get_digraph_scc(diff_g);
    conn = get_pairs_connectivity(scc);
    min_conn = conn;
	if (max_set != 0){
    	memcpy(opt_set[crtl_set_idx], set_gen, k*sizeof(vertex_type));
	}
    n_crtl_set++;
    crtl_set_idx++;
    free_digraph(diff_g);
    free_cc(scc);
	idx_s++;
	
	if ( k != n ){ 
        while( (set_gen[l]+(k-l) <= n || l != 0)){ /* while there are valid elments to increment */
            set_gen[lst_e]++;
            if( set_gen[lst_e] > n ){ /* increment the next incremental element */
                for ( j = lst_e; j > l; j--) { 
                    /* Check if there is an incrementable element between <l> and the last element
                     * of <set_gen> starting from the end of <set_gen>, otherwise we use <l> */
                    if ( set_gen[j]+(k-j) <= n ){
                        l = j;
                        break;
                    }
                }
                if ( set_gen[l]+(k-l) > n ){ /* Check if element <l> can be incremented, otherwise we pass to the next <l--> */
                    l--;
                }
                set_gen[l]++;
                for ( j = l+1; j < k; j++){ /* set the values of the new starting indexes from <l> */
                    set_gen[j] = set_gen[j-1]+1;
                }	
            }
            diff_g = delete_nodes_from_digraph( g, set_gen, k );
            scc = get_digraph_scc(diff_g);
            conn = get_pairs_connectivity(scc);
            if(conn < min_conn){ /* Save the set as the new optimal solution */
                min_conn = conn;
                n_crtl_set = 1;
                crtl_set_idx = 0;
                if (max_set != 0){
                    memcpy(opt_set[crtl_set_idx], set_gen, k*sizeof(vertex_type));
                    crtl_set_idx++;
                }
            }else if ( min_conn == conn ){ /* Add the set to the optimal solution */
                n_crtl_set++;
                if (crtl_set_idx < max_set){
                    memcpy(opt_set[crtl_set_idx], set_gen, k*sizeof(vertex_type));
                    crtl_set_idx++;
                }
            }
            free_digraph(diff_g);
            free_cc(scc);
            idx_s++;
        }
        free(set_gen);
    }
	*end_conn = min_conn;
	*tot_sets = n_crtl_set;

	return opt_set;
}

