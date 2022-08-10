//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

/* File: graph_functions.c */

#include "graph_functions.h"
#include "timer.h"

#define BUFFSIZE     1024
#define ALLOC_BLOCK (2*BUFFSIZE)
/*
 * Function: 	read_graph_rmat
 * --------------------------------------
 * Reads the number of vertices and edges (first line in file) 
 * and the edge list from the input file;
 * it stores "nvertices", "nedges" and return the array "edges"
 *
 * @param fname		the input filename
 * @param nvertices	it will contain the number of vertices
 * @param nedges	it will contain the number of edges
 * @param read_edges	if 0 it will only read nvertices and nedges 
 *
 * @return		the edges list, the number of vertices and edges
 *
 */ 
static VTYPE *read_graph_rmat(const char *fname, VTYPE *nvertices, VTYPE *nedges, int read_edges) {
	FILE *fp;
	char buffer[BUFFSIZE];
	VTYPE buflen;
	VTYPE comment_count;
	VTYPE nv, ne, i, j;
	VTYPE line_count, n;
	VTYPE nmax;
	VTYPE *ed;
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
	ed = (VTYPE *)Malloc(nmax * sizeof(VTYPE));


	while (1) {

		// READ LINES
		fgets(buffer, BUFFSIZE, fp);
		buflen = strlen(buffer);
		if (buflen >= BUFFSIZE) {
			fprintf(stderr, "The line is to long, increase the BUFFSIZE! Exit\n");
			exit(EXIT_FAILURE);
		}

		if (feof(fp)) {
			break;
		} else if (ferror(fp)) {
			fprintf(stderr, "\nAn error ocurred while reading the file\n");
			perror("MAIN:");
		}

		// SCAN THE LINE
		if (strstr(buffer, "#") != NULL) {
			if (strstr(buffer, "Nodes:")) {
				sscanf(buffer, "# Nodes: "TYPE_FORMAT" Edges: "TYPE_FORMAT"\n", &nv, &ne);
				if(!read_edges){
					*nvertices=nv;
					*nedges=ne;
					Free(ed);
					fclose(fp);
					return NULL;
				}
			}
			comment_count++;
		} else {
			if ((nv == 0) || (ne == 0)) {
				fprintf(stderr, "Error reading the number of vertices or edges in %s\n", fname);
				fprintf(stderr, "nv = "TYPE_FORMAT" ne = "TYPE_FORMAT"\n", nv, ne);
				exit(EXIT_FAILURE);
			}

			// Read edges
			if (sscanf(buffer, ""TYPE_FORMAT" "TYPE_FORMAT"\n", &i, &j) != 2) {
				fprintf(stderr, "Error reading graph arc (%s).\n", fname);
				exit(EXIT_FAILURE);
			}

			if (i >= (nv + 1) || j >= (nv + 1)) {
				fprintf(stderr, "found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
						fname, (comment_count + line_count), i, j);
				exit(EXIT_FAILURE);
			} else if ((i == 0) || (j == 0)) {
				fprintf(stderr, "found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
						fname, (comment_count + line_count), i, j);
				fprintf(stderr, "Label 0 is not allowed, labels must start at 1.");
				exit(EXIT_FAILURE);
			}

			if (n >= nmax) {
				nmax += ALLOC_BLOCK;
				ed   = (VTYPE *)Realloc(ed, nmax * sizeof(VTYPE));
			}
			if ( i != j) {
				ed[n]   = i;
				ed[n + 1] = j;
				n       += 2;
			} else {
				self_loops++;
			}
		}
		line_count++;

	}



	fclose(fp);
	if (ne != ((n / 2) + self_loops)) {
		fprintf(stderr, "Error reading the input file %s: the number of edges read differ from the number of edges VTYPE the header\n", fname);
		fprintf(stderr, "nedges header = "TYPE_FORMAT" edges lines read = "TYPE_FORMAT"\n", ne, n / 2);
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


/*
 * Function: read_graph_dimacs
 * -----------------------------------------
 * It reads flowgraph file type DIMACSGRAPH, source at first line of file followed by the edge list,
 *
 * @param fname		the input filename
 * @param nv		it will contain the number of vertices
 * @param ne 		it will contain the number of edges
 * @param root		the root vertex
 * @param read_edges	if 0 it will only read nvertices and nedges 
 *
 * @return 		the edges list, the number of vertices and edges
 *
 */
static VTYPE *read_graph_dimacs(const char *fname, VTYPE *nvertices, VTYPE *nedges, VTYPE *root, int read_edges) {
	VTYPE  *ledges;                       // The edge list read from file
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
			if (sscanf(line, "p "TYPE_FORMAT" "TYPE_FORMAT" "TYPE_FORMAT"", nvertices, nedges, root) != 3) {
				fprintf(stderr, "Error reading graph size (%s).\n", fname);
				exit(EXIT_FAILURE);
			}
			if(!read_edges){
				return ledges;
			}
			// Allocate memory for the edge list
			ledges = Malloc(2 * (*nedges) * sizeof(ledges));

			++nlines;
			continue;
		}


		if (sscanf(line, "a "TYPE_FORMAT" "TYPE_FORMAT"", &u, &v) != 2) {
			fprintf(stderr, "Error reading graph arc (%s).\n", fname);
			exit(EXIT_FAILURE);
		}

		if (u >= (*nvertices + 1) || v >= (*nvertices + 1)) {
			fprintf(stderr, "found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
					fname, nlines, u, v);
			exit(EXIT_FAILURE);
		} else if ((u == 0) || (v == 0)) {
			fprintf(stderr, "found invalid edge in %s, line "TYPE_FORMAT", edge: ("TYPE_FORMAT", "TYPE_FORMAT")\n",
					fname, nlines, u, v);
			fprintf(stderr, "Label 0 is not allowed, labels must start at 1.");
			exit(EXIT_FAILURE);
		}

		if (u == v) {
			self_loops++;
		} else {
			ledges[p++] = u;
			ledges[p++] = v;
			++nlines;
		}

		if (p > (2 * (*nedges))) {
			fprintf(stderr, "Error! Graph has > "TYPE_FORMAT" arcs.\n", (*nedges));
			exit(-1);
		}
	}

	if (self_loops > 0) {
		fprintf(stderr, "[WARNING]: there were "TYPE_FORMAT" self loops\n", self_loops);
		*nedges = p/2;
	}

	fclose(fin);

	return ledges;
}
#undef BUFFSIZE

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
VTYPE *read_graph(int type, const char *fname, VTYPE *n, VTYPE *m, VTYPE *root, int read_edges) {

	VTYPE *edges = NULL;

	switch (type) {

		case DIMACS:
			edges = read_graph_dimacs(fname, n, m, root, read_edges);
			break;

		case RMAT:
			edges = read_graph_rmat(fname, n, m, read_edges);
			break;

		default:
			fprintf(stderr, "Error, UNKNOWN input graph type %d\n", type);
			exit(EXIT_FAILURE);
			break;

	}

	return edges;
}

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
GRAPH *build_graph_datastruct(VTYPE *edges, VTYPE n, VTYPE m, VTYPE root) {
	GRAPH *g;
	VTYPE u, v;
	VTYPE *gin, *g1in, *gout, *g1out;
	VTYPE *g_next_in, *g_next_out;
	VTYPE i;
	VTYPE *labels;
	VTYPE *removed;

	// Allocate memory for the adjacency list
	g     = Calloc(1,      sizeof(GRAPH));
	g1in  = Malloc((n + 2) * sizeof(VTYPE));
	g1out = Malloc((n + 2) * sizeof(VTYPE));
	gin   = Malloc(    m * sizeof(VTYPE));
	gout  = Malloc(    m * sizeof(VTYPE));
	labels = Malloc((n + 1) * sizeof(VTYPE));
	removed = NULL;
	// Auxiliary array
	g_next_in  = Malloc((n + 2) * sizeof(VTYPE));
	g_next_out = Malloc((n + 2) * sizeof(VTYPE));

	memset(g, 	       0,         sizeof(GRAPH));
	memset(g1in,       0, (n + 2)*sizeof(VTYPE));
	memset(g1out,      0, (n + 2)*sizeof(VTYPE));
	memset(g_next_in,  0, (n + 2)*sizeof(VTYPE));
	memset(g_next_out, 0, (n + 2)*sizeof(VTYPE));


	// Read vertices from edge list and store them in gin and gout
	for (i = 0; i < m; i++) {
		u = edges[2 * i];
		v = edges[2 * i + 1];
		g1out[u + 1]++;
		g1in[v + 1]++;
	}

	for (i = 1; i <= (n + 1); i++) {
		g1out[i]      += g1out[i - 1];
		g_next_out[i]  = g1out[i];
		g1in[i]       += g1in[i - 1];
		g_next_in[i]   = g1in[i];
	}

	for ( i = 0; i <= n; i++) {
		labels[i] = i;
	}
	for (i = 0; i < m; i++) {
		u = edges[2 * i];
		v = edges[2 * i + 1];
		gout[g_next_out[u]++] = v;
		gin[g_next_in[v]++]   = u;
	}

	g->gin   = gin;
	g->gout  = gout;
	g->g1in  = g1in;
	g->g1out = g1out;
	g->labels = labels;
	g->n     = n;
	g->m     = m;
	g->root  = root;
	g->removed = removed;
	g->removed_num = 0;
	g->original_n = n;
	g->cn_equiv_len = 0;
	g->cn_equiv = NULL;

	Free(g_next_in);
	Free(g_next_out);
	return g;
}


/*
 * Function: get_reverse_graph
 * -----------------------------------
 * It returns a pointer to the reverse graph
 *
 * @param g	the input graph
 *
 * @return	the reverse graph of g
 *
 */
static GRAPH *get_reverse_graph(GRAPH *g) {
	GRAPH *gR;

	gR        = Malloc(sizeof(GRAPH));

	gR->g1in  = g->g1out;
	gR->g1out = g->g1in;
	gR->gin   = g->gout;
	gR->gout  = g->gin;
	gR->labels = g->labels;
	gR->n = g->n;
	gR->m = g->m;
	gR->root = g->root;
	gR->removed = g->removed;
	gR->removed_num = g->removed_num;
	gR->original_n = g->original_n;
	gR->cn_equiv = g->cn_equiv;
	gR->cn_equiv_len = g->cn_equiv_len;


	return gR;
}


/*
 * Function init_data_structures
 * --------------------------------------
 * It initializes the arrays, labels starts from 1 so when use N = n+1
 *
 * @param tree		the Deep-First-Search tree datastructure	
 * @param clists	The C-lists datastructures
 * @param n 		the number of vertices
 * @param m		the number of edges
 *
 * @return		it modifices the input parameters tree and clisss
 *
 */
static void init_data_structures(DFS_TREE **tree, C_EDGE_LISTS **clists, VTYPE n, VTYPE m) {
	VTYPE mnt = m - n + 2;
	VTYPE N   = n + 1;

	DFS_TREE *t;
	t = Malloc(sizeof(DFS_TREE));
	t->dfsl = Malloc(N * sizeof(VTYPE));
	t->ptol = Malloc(N * sizeof(VTYPE));
	t->p    = Malloc(N * sizeof(VTYPE));
	memset(t->dfsl, 0, N * sizeof(VTYPE));
	memset(t->ptol, 0, N * sizeof(VTYPE));
	memset(t->p,    0, N * sizeof(VTYPE));


	C_EDGE_LISTS *c;
	c = Malloc(sizeof(C_EDGE_LISTS));
	c->cycle = Malloc(mnt * sizeof(VTYPE));
	c->cross = Malloc(mnt * sizeof(VTYPE));
	c->next  = Malloc(mnt * sizeof(VTYPE));
	c->last_cycle = Malloc(N * sizeof(VTYPE));
	c->last_cross = Malloc(N * sizeof(VTYPE));
	memset(c->cycle, 0, mnt * sizeof(VTYPE));
	memset(c->cross, 0, mnt * sizeof(VTYPE));
	memset(c->next,  0, mnt * sizeof(VTYPE));
	memset(c->last_cycle,  0, N * sizeof(VTYPE));
	memset(c->last_cross,  0, N * sizeof(VTYPE));

	c->lastpos = 1;

	(*tree)   = t;
	(*clists) = c;
}

/*
 * Function: push_clists
 * ---------------------------
 * It pushes the vertex s, parent of v into clist of v.
 * clists is organised in three arrays and one integer:
 * 1) array cl->cycle is the array wich actually contains the lists of vertices
 * 2) array cl->last_cycle is an array that in position v contains the pointer
 * to the beginning of the list of vertices in cycle that "belongs" to v
 * 3) array cl->next contatains the indices of the next position of current position
 * 4) value lastpos store the last (current) accessed position in cl->cycle
 * s is the source of the edges, the vertex from wich the edge exit.
 *
 * @param v	
 * @param s
 * @param cl
 *
 * @return 
 *
 */
static void push_clists(VTYPE v, VTYPE s, C_EDGE_LISTS *cl) {
	VTYPE oldp;                     // old position
	VTYPE curp = cl->lastpos;      // current position is stored in the clists data structure - IT CAN BE UNITIALIZED!

	// If C-list is empty, push first element
	if (!cl->last_cycle[v]) {
		cl->cycle[curp]    = s;
		cl->next [curp]    = curp;
		cl->last_cycle[v]  = curp;
	} else {
		oldp = cl->last_cycle[v];
		cl->cycle[curp]   = s;
		cl->next [curp]   = cl->next[oldp];
		cl->next [oldp]   = curp;
		cl->last_cycle[v] = curp;
	}
	curp++;
	cl->lastpos = curp;
}


/*
 * Function: push_cross
 * -----------------------
 * It pushes cross edges
 * Here we write different elements of the cycle arrays but we can use
 * it because the total number of cross + cycle is < m - n + 1
 * Moreover later we will merge lists which is easier within the same array
 * Work a the code above but here I need to push the node t in the list of
 * lca which is the least common ancestor computed via Binary Search
 *
 * @param lca
 * @param v
 * @param s
 * @param cl
 *
 * @return
 *
 */
static void push_cross(VTYPE lca, VTYPE v, VTYPE s, C_EDGE_LISTS *cl) {

	VTYPE curp = cl->lastpos;
	if ( !cl->last_cross[lca] ) {                    // C-list is empty, push first element
		cl->cycle[curp]      = s;
		cl->cross[curp]      = v;
		cl->next[curp]       = curp;
		cl->last_cross[lca]    = curp;
	} else {
		int old_pos = cl->last_cross[lca];
		cl->cycle[curp]     = s;
		cl->cross[curp]     = v;
		cl->next[curp]      = cl->next[old_pos];
		cl->next[old_pos]   = curp;
		cl->last_cross[lca]   = curp;
	}
	curp++;
	cl->lastpos = curp;

}



/*
 * Function: binary_search
 * --------------------------------------
 * It implements the binary search 
 *
 * @param v		the array of num values
 * @param num		the len of v
 * @param val		the value to search
 *
 * @return		the index of the searched value 
 *
 */
static VTYPE binary_search(VTYPE *v, VTYPE num, VTYPE val) {

	VTYPE  min = 0;
	VTYPE  max = num - 1;
	VTYPE  mid = max >> 1;

	while (min <= max) {
		if (v[mid]  < val)      min = mid + 1;
		else                    max = mid - 1;
		mid = (max >> 1) + (min >> 1) + ((min & max) & 1); //(max + min) >> 1
	}
	return mid;
}

/*
 * Function: check_degree
 * -------------------------------------
 * It checks if both the incident and outgoing vertex degree > 0
 * This function is used to check if the choosed root vertex has degree > 0
 *
 * @param v		the node to check
 * @param g		the input graph
 *
 * @return		it returns -1 if the in/out degree of v is 0, otherwise returns 0
 *
 */
static int check_degree(VTYPE v, GRAPH *g) {
	VTYPE outdeg, indeg;
	indeg  = g->g1in[v + 1]  - g->g1in[v];
	outdeg = g->g1out[v + 1] - g->g1out[v];

	if (indeg <= 0 || outdeg <= 0) {
		//    fprintf(stderr, "check_degree: degree of vertex: ("TYPE_FORMAT", "TYPE_FORMAT")\n", outdeg, indeg);
		return -1;
	} else {
#ifdef DEBUG
		printf("vertex degree, out: "TYPE_FORMAT", in :"TYPE_FORMAT"\n", outdeg, indeg);
#endif
	}

	return 0;
}

/*
 * Function: idfs_for_lnt
 * ----------------------------------------------------
 * It computes the DFS plus the list of incoming edges ordered with respect to
 * LeastCommonAncestor. This is used later in the computation of headers
 * It returns the last label (max label) of the DFS tree
 * Note: the variables are unsigned that starts from 1 and 0 is reserved for
 * uninitialised. This is possible because graph labels starts from 1.
 *
 * @param root		the root vertex
 * @param g		the graph
 * @param tree		the DFS tree w.r.t. to g
 * @param clists	the C-lists
 *
 * @return		the last label (max label) of the DFS tree
 *
 */
static VTYPE idfs_for_lnt(VTYPE root, GRAPH *g, DFS_TREE *tree, C_EDGE_LISTS *clists) {
	VTYPE sid, opvid;                // stack last index, index of last opened vertex
	VTYPE i, k;
	VTYPE u, v, w;                   // vertices
	VTYPE ulabel, vlabel;            // vertices labels in dfs

	VTYPE stacksize;
	VTYPE *stack, *lastdesc, *openv;    // stack, array of last descendant of u, array of open vertices

	VTYPE m = g->m;
	VTYPE n = g->n;

	CHECK_GRAPH(g, "idfs_for_lnt");
	int root_degree = check_degree(root, g);
	if ( root_degree < 0 ) {
		FATAL_ERROR("root degree must be greater than 0\n");
	}
	stacksize = 2 * (m + 1) + 2; //very worst case 2*(M+1)!!!
	stack     = Malloc(stacksize * sizeof(VTYPE));
	lastdesc  = Malloc((n + 1)   * sizeof(VTYPE));
	openv     = Malloc(stacksize * sizeof(VTYPE));


	memset(stack, 0, stacksize * sizeof(VTYPE));
	memset(openv, 0, stacksize * sizeof(VTYPE));
	memset(lastdesc, 0, (n + 1)*sizeof(VTYPE));

	// We use a stack to store the current open vertices.
	// To optimise memory accesses the vertex and its predecessor
	// are stored in adjacent locations of the stack (2i, 2i+1).
	// Moreover we store directly the dfs labels of the predecessors.
	// Initialise the stack with the root, dfsl(root) == p(root) = 1
	stack[2] = root;     // sid starts from 1 thus stack starts from 2
	stack[3] = 1;
	sid      = 1;
	opvid    = 0;
	openv[0] = 0;

	// DFS LABELS START FROM 1
	ulabel    = 0;


	while (sid > 0) {
		u = stack[2 * sid];
		w = stack[2 * sid + 1];
		// If u is not yet visited
		if (tree->dfsl[u] == 0) {
			// Labels starts from 1
			ulabel++;
			opvid++;
			// set the label of u and p(u)
			tree->dfsl[u]      = ulabel;
			tree->ptol[ulabel] = u;
			tree->p[ulabel]    = w;       // Here we store the predecessors dfs label (which we need in the headers computation)
			openv[opvid]       = ulabel;  // The last open vertex is inserted in openv


			// Visit the adjacency list of u from the last entry (this way it's
			// easier to preserve the preorder of DFS)
			for (i = (g->g1out[u + 1] - 1); i >= g->g1out[u]; i--) {
				v = g->gout[i];
				vlabel = tree->dfsl[v];
				// If it's the first time we see v we add it to the stack
				if (vlabel == 0) {
					sid++;
					CHECK_STACK_OVERFLOW
						stack[2 * sid]   = v;
					stack[2 * sid + 1] = ulabel;
					// If v was already seen
				} else {
					// Push v label in the Clist
					if (lastdesc[vlabel] == 0) {
						push_clists(vlabel, ulabel, clists);
					} else {
						k = binary_search(openv, opvid, vlabel);
						push_cross(openv[k], vlabel, ulabel, clists);
					}
				}
				// i can be unsigned!!!
				if (i == 0) break;
			}
			// u it's visited
		} else {
			//We need the dfs_label of the current vertex, otherwise we update every time the last labeled vertex
			if (lastdesc[tree->dfsl[u]] == 0) {
				lastdesc[tree->dfsl[u]] = ulabel;
				opvid--;                     // delete vertex in openv
			}
			sid--;                      // delete vertex in stack
		}
	}

	//Memory cleanup
	Free(openv);
	Free(lastdesc);
	Free(stack);

	return ulabel;
}

/*
 * Function: compress
 * --------------------------
 *  Auxliary function used in the process of finding dominators
 *
 *  @param v
 *  @param lparent
 *  @param semi
 *  @param label
 *  @param c
 *  @param stackz
 *
 *  @return 
 *
 */
static VTYPE compress(VTYPE v, VTYPE *lparent, VTYPE *semi, VTYPE *label, VTYPE c, VTYPE *stackz) {
	VTYPE p;
	int  i, j;

	i = 0;
	while (((p = lparent[v]) > c) && (p != v)) {
		stackz[2 * i]   = v;
		stackz[2 * i + 1] = p;
		v = p;
		i++;
	}

	i--;
	for (j = i; j >= 0; j--) {
		v = stackz[2 * j];
		p = stackz[2 * j + 1];
		if (semi[label[p]] < semi[label[v]]) {
			label[v] = label[p];
		}
		lparent[v] = lparent[p];
	}
	return label[v];
}


/*
 * Function: build_dom
 * ----------------------------------
 * It implements the iterative simple Lengauer Tarjan method to find dominators
 *
 * @param g		the input graph
 * @param tree		the DFS tree
 * @param dom_tree	the dominator tree, filled by the function
 *
 * @return		it returns 0, success
 *
 */
static int build_dom(GRAPH *g, DFS_TREE *tree, VTYPE *dom_tree) {
	VTYPE i, j;
	VTYPE u, v, w, s;
	VTYPE n = g->n, m = g->m;
	VTYPE stacksize = 2 * (m + 1) + 2; //very worst case 2*(M+1)!!!

	VTYPE *	semi 		= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE)) ;
	VTYPE * bucket 		= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE));
	VTYPE * label 		= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE));
	VTYPE * ancestor 	= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE));
	VTYPE * stack 		= (VTYPE*)Malloc(stacksize * sizeof(VTYPE));

	for ( i = g->n; i >= 1; i--) {
		semi[i] = label[i] = i;
		ancestor[i] = tree->p[i];
	}
	memset(bucket,	0, (n + 2)*sizeof(VTYPE));
	memset(stack, 	0, stacksize * sizeof(VTYPE));
	// process the vertex in reverse order w.r.t. DFS preorder
	for (i = g->n; i > 1 ; i--) {
		for (v = bucket[i]; v ; v = bucket[v]) {
			u = compress(v, ancestor, semi, label, i, stack);
			dom_tree[v] = (semi[u] < semi[v]) ? u : i;
		}
		// process all the incoming edges of w
		w = tree->ptol[i];
		for (j = g->g1in[w]; j < g->g1in[w + 1]; j++) {
			v = tree->dfsl[g->gin[j]];
			u = (v <= i) ? v : (compress(v, ancestor, semi, label, i, stack));
			if (semi[u] < semi[i]) {
				semi[i] = semi[u];
			}
		}
		s = semi[i];
		if (s != ancestor[i]) { //if semidominator n not lparent: add i to s's bucket
			bucket[i] = bucket[s];
			bucket[s] = i;
		} else {
			dom_tree[i] = s; //semidominator is lparent: s is a candidate dominator
		}
	}
	// process bucket 1
	for (v = bucket[1]; v; v = bucket[v]) {
		dom_tree[v] = 1;
	}

	// recover idoms
	dom_tree[g->root] = 0;

	// make relative absolute
	for (i = 2; i <= g->n; i++) {
		if (dom_tree[i] != semi[i]) dom_tree[i] = dom_tree[dom_tree[i]];
	}
	Free(semi);
	Free(bucket);
	Free(label);
	Free(ancestor);
	Free(stack);
	return 0;
}

/*
 * Function: C_pop
 * ----------------------------------
 * It gets the first vertex in the list of C
 * 
 * @param node
 * @param clists
 *
 * @return		the first vertex inthe list of C
 *
 */
static VTYPE C_pop(VTYPE node, C_EDGE_LISTS *clists) {
	VTYPE source;
	// C-list contains one element; pop last element
	if (clists->next[clists->last_cycle[node]] == clists->last_cycle[node]) {
		source = clists->cycle[clists->last_cycle[node]];
		clists->last_cycle[node] = 0;
	} else {
		source = clists->cycle[clists->next[clists->last_cycle[node]]];      // first element in C-list
		clists->next[clists->last_cycle[node]] = clists->next[clists->next[clists->last_cycle[node]]];
	}

	return source;
}

/*
 * Function: C_popCrooss_tarpos
 * ---------------------------------------------
 * It pops the element off C-list(node) save element and element position in clists->cycle
 *
 * @param node
 * @param arcvertices
 * @param clists
 *
 * @return 
 *
 */
static VTYPE C_popCross_tarpos(VTYPE node, VTYPE *arcvertices, C_EDGE_LISTS *clists) {
	VTYPE source_vertex;
	VTYPE target_vertex;
	VTYPE position;
	if ( clists->next[clists->last_cross[node]] == clists->last_cross[node] ) {   // C-list contains one element
		source_vertex = clists->cycle[clists->last_cross[node]];
		target_vertex = clists->cross[clists->last_cross[node]];
		position = clists->last_cross[node];
		clists->last_cross[node] = 0;
	} else {
		source_vertex = clists->cycle[clists->next[clists->last_cross[node]]];
		target_vertex = clists->cross[clists->next[clists->last_cross[node]]];
		position = clists->next[clists->last_cross[node]];
		clists->next[clists->last_cross[node]] = clists->next[clists->next[clists->last_cross[node]]];
	}
	arcvertices[0] = source_vertex;
	arcvertices[1] = target_vertex;
	return position;
}

/*
 * Function: C_merge
 * -----------------------------
 * It merges the clists of c_n1 and c_n2
 *
 * @param c_n1
 * @param c_n2
 * @param clists
 *
 * @return 
 *
 */
static void C_merge(VTYPE c_n1, VTYPE c_n2, C_EDGE_LISTS *clists) {

	VTYPE first1;
	VTYPE first2;

	if (!clists->last_cycle[c_n2]) return;

	if (!clists->last_cycle[c_n1]) {
		clists->last_cycle[c_n1] = clists->last_cycle[c_n2];
	} else {
		// add second list at the beginning (befor first node) of first
		first1 = clists->next[clists->last_cycle[c_n1]];
		first2 = clists->next[clists->last_cycle[c_n2]];
		clists->next[clists->last_cycle[c_n1]] = first2;
		clists->next[clists->last_cycle[c_n2]] = first1;
	}
	clists->last_cycle[c_n2] = 0;
}

/*
 * Function: C_move2cycle
 * -------------------------
 * It pushes the element in clists->cycle[target_pos] to C-list(Cnode); previously popped by C_pop_tarpos
 *
 * @param node
 * @param target_pos
 * @param clists
 *
 * @return
 *
 */
static void C_move2cycle(VTYPE node, VTYPE target_pos, C_EDGE_LISTS *clists) {
	if ( !clists->last_cycle[node] ) {               // C-list is empty; clists->cycle[tar_pos] first element
		clists->next[target_pos]   = target_pos;
		clists->last_cycle[node]       = target_pos;
	} else {
		int old_pos = clists->last_cycle[node];

		clists->next[target_pos]   = clists->next[old_pos];
		clists->next[old_pos]      = target_pos;
		clists->last_cycle[node]       = target_pos;
	}
}

/*
 * Function find
 * ---------------
 * It returns the header of the biggest loop that contains p
 *
 * @param p
 * @param ufparent
 * @param z
 *
 * @return	the header of the biggest loop that contains p
 *
 */
static VTYPE find(VTYPE p, VTYPE *ufparent, VTYPE *z) {
	VTYPE k = 0;
	VTYPE i;

	while (ufparent[p] != p) {
		z[k] = p;
		p    = ufparent[p];
		k++;
	}
	for (i = 0; i < k; i++) {
		ufparent[z[i]] = p;
	}
	return p;
}

/*
 * Function: build_header
 * ---------------------------
 * It returns the header of each loop
 *
 *  @param g		the input graph
 *  @param tree		the DFS tree
 *  @param clists	the clists
 *  @param cscc_r	
 *
 *  @return	
 *
 */
static VTYPE * build_header(GRAPH *g, DFS_TREE *tree, C_EDGE_LISTS *clists, VTYPE **cscc_r) {
	VTYPE i, n;
	VTYPE v, v_in, tar_pos;
	VTYPE uF;
	VTYPE pointer[2];
	VTYPE sdom;

	n = g->n;

	VTYPE * header	 	= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE));
	VTYPE * z		 	= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE));
	VTYPE * ufparent 	= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE)); //union-find
	VTYPE * cscc	 	= (VTYPE*)Malloc((n + 2) * sizeof(VTYPE)); //cardinality of loop headed by each vertex
	sdom = 0;
	memset(header,	0, (n + 2)*sizeof(VTYPE));
	memset(z,		0, (n + 2)*sizeof(VTYPE));
	for ( i = 0; i < (n + 2); i++) {
		ufparent[i] = i;
		cscc[i] = 1;
	}
	pointer[0] = 0;
	pointer[1] = 0;

	// process the vertices in reverse preorder
	for (i = n; i > 0; i--) {
		sdom = 0;
		while ( clists->last_cross[i] ) {
			// pop from C-list(i); if descendant then unite, if not move to C-cross-list
			tar_pos = C_popCross_tarpos(i, pointer, clists);
			v = pointer[0];	      // source edge
			v_in = pointer[1];        // target edge
			v_in = find(v_in, ufparent, z);
			// move source = clists->cycle[tar_pos] to C_cycle-list(find(v_in)) not yet in loop
			C_move2cycle(v_in, tar_pos, clists);
		}
		// compute vertices added to C_cycle-list[i] while contracting crossedges
		while ( clists->last_cycle[i] ) {
			// pop from C_cycle-list(i)
			v = C_pop(i, clists);
			while ( (uF = find(v, ufparent, z)) != i ) {
				C_merge(i, uF, clists);     // add C-list(uF) of cross-arcs incident in uF to C_cycle-list(i); clists->last_cycle(i)
				header[uF] = i;
				ufparent[uF] = i;
				cscc[i] += cscc[uF];
				v = tree->p[uF];
				sdom++;
			}
		}
	}
	header[g->root] = sdom;
	Free(z);
	Free(ufparent);
	*cscc_r = cscc;
	return header;
}


/*
 * Function: remove_vertex_from_graph
 * --------------------------------------
 *  It removes a selected vertex from a graph
 *
 *  @param g		the input graph
 *  @param v		the vertex to be removed
 *
 *  @return		the graph G/v
 *
 */
static GRAPH * remove_vertex_from_graph(GRAPH *g, VTYPE v) {
	VTYPE i;					/* Index variable */
	VTYPE j;					/* Index variable */
	VTYPE k;					/* Index variable */
	VTYPE h;					/* Index variable */
	VTYPE removed;				/* Number of removed edges in G after deleting v */
	GRAPH *g_minus_v;			/* The new graph G\v */

	removed = (g->g1out[v + 1]) - (g->g1out[v]) + (g->g1in[v + 1]) - (g->g1in[v]);
	g_minus_v = (GRAPH*) Calloc( 1,sizeof(GRAPH));
	g_minus_v->n = g->n - 1 ;
	g_minus_v->m = g->m - removed;
	g_minus_v->root = 1;
	g_minus_v->removed_num = g->removed_num;
	g_minus_v->original_n = g->original_n;

	g_minus_v->labels = Malloc( (g_minus_v->n + 1)*sizeof(VTYPE));
	g_minus_v->g1in = Malloc( (g_minus_v->n + 2)*sizeof(VTYPE));
	g_minus_v->g1out = Malloc( (g_minus_v->n + 2)*sizeof(VTYPE));
	g_minus_v->gin = Malloc( (g_minus_v->m + 1) * sizeof(VTYPE));
	g_minus_v->gout = Malloc( (g_minus_v->m + 1) * sizeof(VTYPE));
	g_minus_v->removed = Malloc( (g_minus_v->removed_num +1 )*sizeof(VTYPE));
	if ( g->removed_num > 0){
		memcpy( g_minus_v->removed,  g->removed, (g->removed_num )*sizeof(VTYPE));
	}
	g_minus_v->removed[g_minus_v->removed_num] = g->labels[v];
	g_minus_v->removed_num++;
	g_minus_v->cn_equiv = NULL;
	g_minus_v->cn_equiv_len = 0;

	for (j = 1, k = 0, h = 1, removed = 0 ; j < (g->n+1); j++) { /* Checks all vertices, scans g1out and gout */
		if ( j != v ) {
			g_minus_v->labels[h] = g->labels[j];
			g_minus_v->g1out[h] = g->g1out[j] - removed;
			h++;
		}
		if (j == v) { /* We remove all edges starting from v */
			removed = removed + g->g1out[j + 1] - g->g1out[j];
		} else {
			VTYPE u;
			for ( i = g->g1out[j]; i < g->g1out[j + 1]; i++) { /* All neighbors <u> of <j> */
				u = g->gout[i];
				if (u == v) { /* <u> is the vertex <v> to remove thus we remove edge (<j>,<v>), i.e. we do not copy <u> in the new graph */
					removed++;
				} else { /* <u> is a vertex we do not remove, we copy it with its new ID id needed */
					if (u > v) {
						g_minus_v->gout[k] = u - 1;
					} else {
						g_minus_v->gout[k] = u;
					}
					k++;
				}
			}
		}
	}/* O(n+m) */
	g_minus_v->g1out[h] = g->g1out[j] - removed; /* Update g1out last index [n+1]*/

	for (j = 1, k = 0, h = 1, removed = 0 ; j < (g->n+1); j++) { /* Checks all vertices, scans g1in and gin */
		if ( j != v ) {
			g_minus_v->g1in[h] = g->g1in[j] - removed;
			h++;
		}
		if (j == v) { /* We remove all edges entering v */
			removed = removed + g->g1in[j + 1] - g->g1in[j];
		} else {
			VTYPE u;
			for ( i = g->g1in[j]; i < g->g1in[j + 1]; i++) { /* All neighbors <u> of <j> */
				u = g->gin[i];
				if (u == v) { /* <u> is the vertex <v> to remove thus we remove edge (<j>,<v>), i.e. we do not copy <u> in the new graph */
					removed++;
				} else { /* <u> is a vertex we do not remove, we copy it with its new ID id needed */
					if (u > v) {
						g_minus_v->gin[k] = u - 1;
					} else {
						g_minus_v->gin[k] = u;
					}
					k++;
				}
			}
		}
	}/* O(n+m) */
	g_minus_v->g1in[h] = g->g1in[j] - removed; /* Update g1out last index [n+1]*/

	return g_minus_v;
}


/*
 * Function: compute_scc_tarjan
 * -------------------------------
 *  It computes and returns the Strongly Connected Components of g
 *
 *  @param g		the input graph
 *
 *  @return		the SCC of g
 *
 */
static SCC * compute_scc_tarjan(GRAPH *g ){ 
	VTYPE *lowlink;					/* Contains the lowlink of each vertex, i.e. the root of the scc to which the vertex belong. */
	VTYPE *number;					/* Contains the dfs number of each vertex */
	VTYPE *scc_stack;				/* Contains the verticies of the current scc */

	unsigned char *vertex_in_scc_stack;		/* Array used to check if a vertex is in the <scc_stack> */

	VTYPE *stack;					/* Contains the verticies to visit */
	VTYPE *frontier;				/* Keeps track of the vertex on top of the branch tree we are visiting */

	VTYPE dfs_number;				/* Preorder number assigned by DFS */
	VTYPE stack_index;				/* Index of the top vertex in the stack */
	VTYPE i;						/* Index variable */
	VTYPE scc_stack_index;			/* Vertices 'possibly' part of the current Strongly Connected Component */
	VTYPE frontier_index;			/* Index of the top vertex in the frontier */
	VTYPE w;						/* The neighbor of v */
	VTYPE v;						/* The vertex we are visting */
	VTYPE n_scc;					/* Number of SCC */

	VTYPE stack_index_start_value;	/* Value of <stack_index> before visiting the neighbors of node <v> */

	VTYPE scc_index;				/* Index used to insert information in the <SCC> data structure */
	SCC *scc;						/* Contains the SCC of the graph <g> */

	VTYPE n;						/* Number of vertecies in <g> */
	VTYPE m;						/* Number of edges in <g> */
	VTYPE last_unnumbered;			/* Last not visited/unnambered vertex */

	n = g->n;
	m = g->m;
	lowlink = Calloc( n+1, sizeof(VTYPE));
	number = Calloc( n+1, sizeof(VTYPE));
	scc_stack = Calloc( n+1, sizeof(VTYPE));
	stack = Calloc( m+n+1, sizeof(VTYPE));	
	frontier = Calloc( n+1, sizeof(VTYPE));
	vertex_in_scc_stack = Calloc( n+1, sizeof(unsigned char));

	scc = (SCC *) Calloc( 1, sizeof(SCC));
	scc->vertex = Calloc( n, sizeof(VTYPE));
	scc->index  = Calloc( n+1, sizeof(VTYPE));

	scc_index = 0;
	n_scc = 0;

	scc_stack_index = 1;
	frontier_index = 1;
	dfs_number = 1;
	stack_index = 1;
	stack[1] = 1;
	last_unnumbered = 1;
	v=0;

	while( (stack_index != 0) || (dfs_number <= n) ){ /* There are still items in the stack or unnumbered (not visited) vertices */
		if ( stack_index != 0 ){ /* There are still vertices to visit in the stack */
			v = stack[stack_index];
		}else{ /* The stack is empty but there are still unnumbered vertices */
			for( i = last_unnumbered; i <= n; i++ ){ /* Select the new vertex for the visit, <stack_index> is 0 */
				if(number[i] == 0){
					v = i;
					stack_index++;
					stack[stack_index] = v;	
					last_unnumbered = i+1;
					break;
				}
			}
		}
		/* <v> is the next vertex to visit (it was in the <stack> or it was still unnumbered) */
		if(number[v] == 0){ /* <v> has not been visited yet */
			lowlink[v] = number[v] = dfs_number;    
			dfs_number++;
			scc_stack[scc_stack_index] = v;			/* we keep track of the vertices in the same component (possibly strongly connected) */
			scc_stack_index++;
			vertex_in_scc_stack[v] = 1;
			frontier[frontier_index] = v;			/* we are visiting <v> and possibly its neighbors */
			frontier_index++;
			stack_index_start_value = stack_index;
			/* Push the neighbors of v in the stack, in reverse order. Only if they are not visited yet. */
			if (g->g1out[v+1] != 0) {
				for ( i = (g->g1out[v+1])-1; i >= g->g1out[v]; i--){
					w = g->gout[i];
					if( number[w] == 0 ){ /* <w> is a neighbor of <v> not yet visited */
						stack_index++;
						stack[stack_index] = w;
					}else{ 
						if(number[w] < number[v]){
							/* We visited <w> before <v> thus v->w is a frond or cross-link.
							 * - Fronds connect descendants with ancestors.
							 * - Cross-Links run from one subtree to another in the tree */
							if ( (vertex_in_scc_stack[w] == 1) 
									&& (number[w] < lowlink[v]) ){ 
								/* check if <w> is in the <scc_stack> (possibly the same scc of <v>) 
								   AND <w> number is lower than the lowlink of <v> */
								lowlink[v] = number[w];
							}
						}
					}
					if (i == 0) {break;} /* Check needed because <i> could be unsigned */
				}
			}
			if( stack_index_start_value == stack_index){ 
				/* <v> is a leaf of the dfs tree, <v> has no unvisited neighbors */
				frontier_index--;
				stack_index--;
				if ( lowlink[v] == number[v] ){ /* Check if <v> is a root of a scc */
					/* The scc of <v> is composed only by <v> */
					n_scc++;
					scc_stack_index--; 
					vertex_in_scc_stack[v] = 0;
					(scc->vertex)[scc_index] = v;
					scc_index++;
					(scc->index)[n_scc] = scc_index;
				}
			}
		}else{/* Vertex <v> has been already visited. */ 
			stack_index--;
			if ( frontier[frontier_index-1] == v ){ /* We are evaluating <v> after the visit of all its descendants, it is the current vertex */
				frontier_index--;
				if (g->g1out[v+1] != 0) { 
					for ( i = (g->g1out[v+1])-1; i >= g->g1out[v]; i--){
						w = g->gout[i];
						if ( (vertex_in_scc_stack[w] == 1) && (lowlink[v] > lowlink[w]) ){ 
							/* <w> is in the <scc_stack> and is a child of <v> in the dfs tree AND <w> lowlink is lower than the lowlink of <v>*/
							lowlink[v] = lowlink[w];
						}
						if (i == 0) {break;} /* Check needed because <i> could be unsigned */
					}
				}
				if ( lowlink[v] == number[v] ){ /* Check if <v> is a root of a scc */
					n_scc++;
					for ( i = scc_stack_index-1; i > 0 ; i--){
						w = scc_stack[i];
						if( number[w] >= number[v]){ /* <w> in the <scc_stack> is a vertex in the same SCC of <v>*/
							scc_stack_index--; /* remove <w> form <scc_stack> */
							vertex_in_scc_stack[w] = 0;
							(scc->vertex)[scc_index] = w;
							scc_index++;
						}else{ break; }
					}
					(scc->index)[n_scc] = scc_index;
				}
			}/* <v> was still in the stack but it is not the current vertex (it was not in the frontier) */
		}
	}
	scc->number = n_scc;
	scc->index = (VTYPE *) Realloc(scc->index, sizeof(VTYPE)*(n_scc+1) );

	Free(lowlink);
	Free(number);
	Free(scc_stack);
	Free(stack);	
	Free(frontier);
	Free(vertex_in_scc_stack);

	return scc;
}


/*
 * Function: free_cc
 * ----------------------------
 *  It deallocates the memory belonging to a CC
 *
 *  @param s 		the input CC
 *
 *  @return
 *
 */
static void free_cc(CC *s) {
	if ( NULL == s) {
		return;
	}

	if ( NULL != s->vertex) {
		Free(s->vertex);
	}

	if ( NULL != s->index) {
		Free(s->index);
	}

	Free(s);
}


/*
 * Function: free_graph
 * -------------------------
 * It deallocates the memory belonging to a graph
 *
 * @param g		the input graph
 *
 * @return
 *
 */
static void free_graph(GRAPH * g) {
	if ( NULL == g) {
		return;
	}

	if (NULL != g->gin) {
		Free(g->gin);
	}

	if (NULL != g->g1in) {
		Free(g->g1in);
	}

	if (NULL != g->gout) {
		Free(g->gout);
	}

	if (NULL != g->g1out) {
		Free(g->g1out);
	}

	if (NULL != g->labels) {
		Free(g->labels);
	}

	if (NULL != g->removed){
		Free(g->removed);
	}
	if (NULL != g->cn_equiv){
		Free(g->cn_equiv);
	}

	Free(g);
}


/*
 * Function: get_art_points
 * ---------------------------------------
 * It return the set of articolation points as the union of the
 * vertices in D(s) and D^R(s)
 * Ref: "Finding strong bridges and strong articulation points in linear time"
 * G. F. Italiano, L. Laura, F. Santaroni
 *
 * @param g		the input graph
 * @param dom_tree	the dominator tree
 * @param dom_treeR	the dominator tree of gR, the reverse graph of g
 * @param num
 *
 * @return		the list of articulation points of g
 *
 */
static VTYPE *get_art_points(GRAPH *g, VTYPE *dom_tree, VTYPE *dom_treeR, VTYPE *num) {
	VTYPE *art_points, *art_points_tmp;
	VTYPE i, j, k;
	VTYPE n = (g->n) + 1;
	VTYPE include_root;

	/*Check root*/ 
	GRAPH* g_minus_root = remove_vertex_from_graph(g,g->root);
	SCC *scc = compute_scc_tarjan(g_minus_root);
	if( scc->number != 1){
		include_root = 1;
	}
	else{
		include_root = 0;
	}
	free_graph(g_minus_root);
	free_cc(scc);

	art_points_tmp = Malloc( n * sizeof(VTYPE));
	memset(art_points_tmp, 0, n * sizeof(VTYPE));

	for ( i = 1, j = 0; i < n; i++ ) {
		if(dom_tree[i] >0 && art_points_tmp[dom_tree[i]] != 1){
			art_points_tmp[dom_tree[i]] = 1;	
			j++;
		}
		if(dom_treeR[i] > 0 && art_points_tmp[dom_treeR[i]] != 1){
			art_points_tmp[dom_treeR[i]]=1;
			j++;
		}
	}
	art_points = Malloc( j * sizeof(VTYPE));
	for(i = 0, k = 0; i < n && k < j; i++){
		if(art_points_tmp[i] != 0){
			if(i != g->root || (i == g->root && include_root) ){
				art_points[k]=i;
				k++;
			}
		}
	}
	Free(art_points_tmp);
	*num = k;
	return art_points;
}


/*
 * Function: translate_dom
 * --------------------------
 * It returns the original node labels for each dominator
 *
 * @param g		the input graph
 * @param tree		the DFS tree
 * @param root		the root vertex
 * @param dom		the list of dominators
 *
 * @return		the original label of each dominator
 *
 */
static VTYPE *translate_dom(GRAPH *g, DFS_TREE *tree, VTYPE root, VTYPE *dom) {
	VTYPE i, n;
	VTYPE *idom = NULL;

	n = (g->n) + 1;
	idom = Malloc(n * sizeof(VTYPE));
	memset(idom, 0, n * sizeof(VTYPE));
	for (i = 2; i < n; i++) {
		idom[tree->ptol[i]] = tree->ptol[dom[i]];
	}
	return idom;
}

/*
 * Function: translate_header
 * --------------------------
 * It returns the original node labels for each header
 *
 * @param g		the input graph
 * @param tree		the DFS tree
 * @param header	the list of headers
 *
 * @return		the original label of each header
 *
 */
static VTYPE *translate_header(GRAPH *g, DFS_TREE *tree, VTYPE *header) {
	VTYPE i, n;
	VTYPE *iheader;

	n = g->n + 1;
	iheader = Malloc(n * sizeof(VTYPE));
	for ( i = 1; i < n; i++) {
		iheader[tree->ptol[i]] = tree->ptol[header[i]];
	}
	return iheader;
}


/*
 * Function: translate_cscc
 * --------------------------
 * It returns the original node labels for each cscc
 *
 * @param g		the input graph
 * @param tree		the DFS tree
 * @param cscc		the list of cscc
 *
 * @return		the original label of each cscc
 *
 */
static VTYPE *translate_cscc(GRAPH *g, DFS_TREE *tree, VTYPE *cscc) {
	VTYPE i, n;
	VTYPE *icscc;

	n = g->n + 1;
	icscc = Malloc(n * sizeof(VTYPE));
	icscc[0] = 0;
	for ( i = 1; i < n; i++) {
		icscc[i] = cscc[tree->dfsl[i]];
	}
	return icscc;
}

/*
 * Function: build_graph_from_array
 * -------------------------------------------------
 * It builds the graph corresponding to Dominator/LoopNesting Tree
 *
 * @param n		the number of vertices
 * @param m		the number of edges
 * @param array		the array containing the dominator/Loop Nesting tree relations
 * @param root		the root vertex
 *
 * @return		the graph object corresponding to the dominator/LoopNesting tree
 *
 */
static GRAPH *build_graph_from_array(VTYPE n, VTYPE m, VTYPE *array, VTYPE root) {
	GRAPH *g;
	VTYPE i;
	VTYPE *gin, *g1in, *gout, *g1out;
	VTYPE *g_next_in, *g_next_out;

	// Allocate memory for the adjacency list
	g     = Calloc(1,          sizeof(GRAPH));
	g1in  = Malloc((n + 2) * sizeof(VTYPE));
	g1out = Malloc((n + 2) * sizeof(VTYPE));
	gin   = Malloc(    m   * sizeof(VTYPE));
	gout  = Malloc(    m   * sizeof(VTYPE));
	// Auxiliary array
	g_next_in  = Malloc((n + 2) * sizeof(VTYPE));
	g_next_out = Malloc((n + 2) * sizeof(VTYPE));

	memset(g1in,       0, (n + 2)  * sizeof(VTYPE));
	memset(g1out,      0, (n + 2)  * sizeof(VTYPE));
	memset(g_next_in,  0, (n + 2)  * sizeof(VTYPE));
	memset(g_next_out, 0, (n + 2)  * sizeof(VTYPE));

	VTYPE u, v;
	for ( i  = 0; i < n ; i++) {
		u = array[i + 1];
		v = i + 1;
		g1out[u + 1]++;
		g1in[v + 1]++;
	}
	for (i = 1; i <= (n + 1); i++) {
		g1out[i]      += g1out[i - 1];
		g_next_out[i]  = g1out[i];
		g1in[i]       += g1in[i - 1];
		g_next_in[i]   = g1in[i];
	}
	for (i = 0; i < n; i++) {
		u = array[i + 1];
		v = i + 1;
		gout[g_next_out[u]++] = v;
		gin[g_next_in[v]++]   = u;
	}
	g->gin   = gin;
	g->gout  = gout;
	g->g1in  = g1in;
	g->g1out = g1out;
	g->n     = n;
	g->m     = m;
	g->labels = NULL;
	g->root  = root;
	g->removed = NULL;
	g->removed_num = 0;
	g->original_n = n;
	g->cn_equiv_len = 0;
	g->cn_equiv = NULL;

	Free(g_next_in);
	Free(g_next_out);
	return g;
}

/*
 * Function: simple_dfs
 * ------------------------------
 * It computes only the DFS
 * It also fills the array of last descendants for each vertex u
 * Note: the variables are unsigned that starts from 1 and 0 is reserved for
 * uninitialised. This is possible because graph labels starts from 1.
 *
 * @param root		the root vertex
 * @param g		the input graph
 * @param lastdesc	the last descendant array
 *
 * @return		the DFS tree
 *
 */
static DFS_TREE *simple_dfs(VTYPE root, GRAPH *g, VTYPE *lastdesc) {
	VTYPE sid;                // stack last index
	VTYPE i;
	VTYPE u, v, w;                   // vertices
	VTYPE ulabel, vlabel;            // vertices labels in dfs

	VTYPE stacksize;
	VTYPE *stack;    // stack

	VTYPE m = g->m;
	VTYPE n = g->n;

	VTYPE N   = n + 1;
	CHECK_GRAPH(g, "simple_dfs");

	int root_degree = check_degree(root, g);
	if ( root_degree < 0) {
		FATAL_ERROR("root degree must be greater than 0\n");
	}

	DFS_TREE *tree;
	tree = Malloc(sizeof(DFS_TREE));
	tree->dfsl = Malloc(N * sizeof(VTYPE));
	tree->ptol = Malloc(N * sizeof(VTYPE));
	tree->p    = Malloc(N * sizeof(VTYPE));
	memset(tree->dfsl, 0, N * sizeof(VTYPE));
	memset(tree->ptol, 0, N * sizeof(VTYPE));
	memset(tree->p,    0, N * sizeof(VTYPE));

	stacksize = 2 * (m + 1) + 2; //very worst case 2*(M+1)!!!
	stack     = Malloc(stacksize * sizeof(VTYPE));


	memset(stack, 0, stacksize * sizeof(VTYPE));
	memset(lastdesc, 0, (n + 1) * sizeof(VTYPE));

	// We use a stack to store the current open vertices.
	// To optimise memory accesses the vertex and its predecessor
	// are stored in adjacent locations of the stack (2i, 2i+1).
	// Moreover we store directly the dfs labels of the predecessors.
	// Initialise the stack with the root, dfsl(root) == p(root) = 1
	stack[2] = root;     // sid starts from 1 thus stack starts from 2
	stack[3] = 1;
	sid      = 1;

	// DFS LABELS START FROM 1
	ulabel    = 0;

	while (sid > 0) {
		u = stack[2 * sid];
		w = stack[2 * sid + 1];
		// If u is not yet visited
		if (tree->dfsl[u] == 0) {
			// Labels starts from 1
			ulabel++;
			// set the label of u and p(u)
			tree->dfsl[u]      = ulabel;
			tree->ptol[ulabel] = u;
			tree->p[ulabel]    = w;       // Here we store the predecessors dfs label (which we need in the headers computation)


			// Visit the adjacency list of u from the last entry (this way it's
			// easier to preserve the preorder of DFS)
			if ( g->g1out[u + 1] != 0) {
				for (i = (g->g1out[u + 1] - 1); i >= g->g1out[u]; i--) {
					v = g->gout[i];
					vlabel = tree->dfsl[v];
					// If it's the first time we see v we add it to the stack
					if (vlabel == 0) {
						sid++;
						CHECK_STACK_OVERFLOW
							stack[2 * sid]   = v;
						stack[2 * sid + 1] = ulabel;
						// If v was already seen
					}
					// i can be unsigned!!!
					if (i == 0) break;
				}
			}
			// u it's visited
		} else {
			if (lastdesc[tree->dfsl[u]] == 0) {
				lastdesc[tree->dfsl[u]] = ulabel;
			}
			sid--;                      // delete vertex in stack
		}
	}

	//Memory cleanup
	Free(stack);

	return tree;
}


/*
 * Function: compute_sumDudNew
 * -----------------------------------------
 *  It computes the total dfo
 *
 *  @param dfo_0
 *  @param dfo_1
 *  @param dfo_01	
 *  @param g 		the input graph
 *  @param w		the array of w
 *  @param out		
 *  @param art_pts	the array of articulation points
 *  @param art_pts_num	the number of articulation points
 *
 *  @return 		the total dfo
 *  
 */
static LTYPE compute_sumDudNew(BUNDLES *dfo_0, BUNDLES *dfo_1, BUNDLES *dfo_01, GRAPH *g, VTYPE *w, LTYPE *out, VTYPE *art_pts, VTYPE art_pts_num) {
	VTYPE i,j, n,deg, tmp_d;
	LTYPE  val, cur;
	n = g->n;
	val = 0;
	deg=0;
	LTYPE sqr = ( (LTYPE)n * (n  - 1)) / 2;
	for ( j = 0 ; j < art_pts_num ; j++) {
		i = art_pts[j];
		VTYPE pca = n - dfo_0->sumC[i] - dfo_1->sumC[i] - 1 + dfo_01->sizeC[i];
		LTYPE pca_sqr = ( (LTYPE)pca * (pca - 1));
		if ( dfo_0->sumC[i] == 0 && dfo_1->sumC[i] == 0)
			continue;
		LTYPE res = 0;
		res += ( (dfo_0->sum_sqC[i] - (dfo_0->sumC[i])) / 2) ;
		res += ( (dfo_1->sum_sqC[i] - (dfo_1->sumC[i])) / 2) ;
		res -=  (dfo_01->sum_sqC[i] - (dfo_01->sumC[i])) / 2;
		if (pca > 1) {
			res += (pca_sqr / 2) ;
		}



		cur = sqr - res;
		if(cur==val){
		    tmp_d=(g->g1out[i+1] - g->g1out[i])	+ (g->g1in[i+1] - g->g1in[i]);
		    if(tmp_d > deg){
			deg=tmp_d;
			*w=i;
		     }
		}else if (cur > val) {
		        deg=(g->g1out[i+1] - g->g1out[i])	+ (g->g1in[i+1] - g->g1in[i]);
			val = cur;
			*w = i;
			*out = res;
		}
		
	}

	return val;
}



#define ELEM_LIST_SIZE 3
/*
 * Function: compute_dfo
 * -------------------------------------
 * f(D~(v)) and f(DR~(v)) for all v in SAP(G)
 * dfs_dom[i] contains the index of the vertex v whose dfs label is i
 *
 * @param g		the input graph
 * @param dom_tree	the dominator tree
 * @param ln_tree	the loop nesting tree
 * @param cscc		the cscc
 * @param dfs_dom	the DFS of dominator tree
 *
 * @return 
 *
 */
static BUNDLES *compute_dfo(GRAPH *g, VTYPE *dom_tree, VTYPE *ln_tree, VTYPE *cscc, VTYPE *dfs_dom) {
	BUNDLES *dfo;
	LTYPE *sumC;
	LTYPE *sum_sqC;
	VTYPE i, n;
	n = g->n;
	dfo 	= Malloc( sizeof(BUNDLES));
	sumC 	= Malloc( (n + 2) * sizeof(LTYPE));
	sum_sqC = Malloc( (n + 2) * sizeof(LTYPE));

	memset(sumC		, 0, (n + 2) * sizeof(LTYPE));
	memset(sum_sqC	, 0, (n + 2) * sizeof(LTYPE));

	for ( i = 2; i <= n ; i++) {
		if (dom_tree[i] != dom_tree[ln_tree[i]]) {
			LTYPE cursq = (LTYPE) cscc[i] * cscc[i];
			sumC[dom_tree[i]] 			 	+= cscc[i];
			sum_sqC[dom_tree[i]] 		 	+=  cursq;
			sumC[dom_tree[ln_tree[i]]] 		-= cscc[i];
			sum_sqC[dom_tree[ln_tree[i]]]  	-= cursq;
		}
	}

	// for each v in D bottom up fashion
	for ( i = n  ; i > 0 ; i--) {
		VTYPE v = dfs_dom[i]; 
		sumC[dom_tree[v]] 		+= sumC[v];
		sum_sqC[dom_tree[v]] 	+= sum_sqC[v];
		// i can be unsigned
		if ( 0 == i ) {
			break;
		}
	}

	dfo->sumC 		= sumC;
	dfo->sum_sqC 	= sum_sqC;
	dfo->sizeC = NULL;
	return dfo;
}


/*
 * Function: build_q
 * ---------------------------------
 * It builds common dominator forest, assign TreeID to each vertex
 * It fills qID and returns Q
 *
 * @param g 		the input graph
 * @param dom		the dominator tree of g
 * @param domR		the dominator tree of gR
 * @param qID		the q ids
 * @param dom_dfs	the DFS of the dominator tree of g
 *
 * @return Q
 *
 */
static VTYPE *build_q(GRAPH *g, VTYPE *dom, VTYPE *domR, VTYPE *qID, DFS_TREE *dom_dfs) {
	VTYPE i, n;
	VTYPE *Q;
	VTYPE treeID;

	n = g->n;
	treeID = 1;

	Q = Malloc((n + 1) * sizeof(VTYPE));
	memset(Q, 0, (n + 1) * sizeof(VTYPE));

	// for each v in D bottom up fashion
	for ( i = n ; i > 1 ; i--) {
		VTYPE v = dom_dfs->ptol[i];
		VTYPE dv = dom[v];
		if (dv == g->root) {
			if (0 == qID[v]) {
				qID[v] = treeID;
				treeID++;
			}
			continue;
		}
		if (v == domR[dv]) {
			Q[v] = dv;
			if ( 0 == qID[v]) {
				qID[v] = qID[dv] = treeID;
				treeID++;
			} else {
				qID[dv] = qID[v];
			}
		} else {
			if ( 0 == qID[v]) {
				qID[v] = treeID;
				treeID++;
			}
			if ( 0 == qID[dv]) {
				qID[dv] = treeID;
				treeID++;
			}

		}
	}
	qID[g->root] = treeID;
	return Q;
}

/*
 * Function: find_z_dfs
 * --------------------------
 * It performs a dfs traversal of dominator tree while it identifies the value of z for every node z where W[v] != null 
 * Let T be a rooted tree and let F be a forest that results from T after deleting some tree edges. 
 * Each tree in F has a unique integer id. We are given a list of pairs (node, TreeID). 
 * For each such pair <x,j> we wish to locate the nearest ancestor z of x in T such that z is in the tree of F with id=j.
 * In our case T=D (the dominator tree), F=Q (the common dominator forest), and for each pair <x,j>, j = TreeID(w(x)).
 *
 * @param g		the input graph
 * @param qID		the q ids
 * @param idom_tree	the dom tree
 * @param root		the root vertex
 *
 * @return
 *
 */
static VTYPE *find_z_dfs(GRAPH *g, VTYPE *W, VTYPE *qID, VTYPE *idom_tree, VTYPE root){
	VTYPE sid;                // stack last index
	VTYPE i;
	VTYPE u, v, w;                   // vertices
	VTYPE ulabel, vlabel;            // vertices labels in dfs

	VTYPE stacksize;
	VTYPE *stack;    // stack
	VTYPE *path;    // path

	VTYPE m;
	VTYPE n;
	VTYPE p_idx;
	VTYPE *z_array = NULL;
	VTYPE *last = NULL;
	VTYPE *tree_id = NULL;


	GRAPH *g_dom = NULL;
	n = g->n;
	m = g->m;
	g_dom = build_graph_from_array(g->n, g->n, idom_tree, root);
	VTYPE N   = n + 1;


	CHECK_GRAPH(g, "simple_dfs");

	int root_degree = check_degree(root, g_dom);
	if ( root_degree < 0) {
		FATAL_ERROR("root degree must be greater than 0\n");
	}

	DFS_TREE *tree;
	tree = Malloc(sizeof(DFS_TREE));
	tree->dfsl = Malloc(N * sizeof(VTYPE));
	tree->ptol = Malloc(N * sizeof(VTYPE));
	tree->p    = Malloc(N * sizeof(VTYPE));
	memset(tree->dfsl, 0, N * sizeof(VTYPE));
	memset(tree->ptol, 0, N * sizeof(VTYPE));
	memset(tree->p,    0, N * sizeof(VTYPE));


	stacksize = 2 * (m + 1) + 2; //very worst case 2*(M+1)!!!
	stack     = Malloc(stacksize * sizeof(VTYPE));
	memset(stack, 0, stacksize * sizeof(VTYPE));


	last       = Malloc(N * sizeof(VTYPE));
	path       = Malloc(N * sizeof(VTYPE));
	tree_id    = Malloc(N * sizeof(VTYPE));
	z_array    = Malloc(N * sizeof(VTYPE));
	memset(last, 0, (n + 1) * sizeof(VTYPE));
	memset(path, 0, (n + 1) * sizeof(VTYPE));
	memset(tree_id, 0, (n + 1) * sizeof(VTYPE));
	memset(z_array,    0, N * sizeof(VTYPE));

	// Filling tree_id
	for( i = 1; i < N; i++){
		// it should be tree_id[v] = qID[W[v]] if W[v] != NULL/0 
		if (W[i] != 0){
			tree_id[i] = qID[W[i]];
		}
	}

	// We use a stack to store the current open vertices.
	// To optimise memory accesses the vertex and its predecessor
	// are stored in adjacent locations of the stack (2i, 2i+1).
	// Moreover we store directly the dfs labels of the predecessors.
	// Initialise the stack with the root, dfsl(root) == p(root) = 1
	stack[2] = root;     // sid starts from 1 thus stack starts from 2
	stack[3] = 1;
	sid      = 1;
	p_idx = 0;

	// DFS LABELS START FROM 1
	ulabel    = 0;

	while (sid > 0) {
		u = stack[2 * sid];
		w = stack[2 * sid + 1];

		// If u is not yet visited
		if (tree->dfsl[u] == 0) {
			// Labels starts from 1
			p_idx++;
			path[p_idx] = u;
			// I must update last for every node (even already visited) 
			last[qID[u]] = u;
			z_array[u] = last[tree_id[u]];

			ulabel++;

			// set the label of u and p(u)
			tree->dfsl[u]      = ulabel;
			tree->ptol[ulabel] = u;
			tree->p[ulabel]    = w;       // Here we store the predecessors dfs label (which we need in the headers computation)


			// Visit the adjacency list of u from the last entry (this way it's
			// easier to preserve the preorder of DFS)
			if ( g_dom->g1out[u + 1] != 0) {
				for (i = (g_dom->g1out[u + 1] - 1); i >= g_dom->g1out[u]; i--) {
					v = g_dom->gout[i];
					vlabel = tree->dfsl[v];
					// If it's the first time we see v we add it to the stack
					if (vlabel == 0) {
						sid++;
						CHECK_STACK_OVERFLOW
							stack[2 * sid]   = v;
						stack[2 * sid + 1] = ulabel;
						// If v was already seen
					}
					// i can be unsigned!!!
					if (i == 0) break;
				}
			}
			// u it's visited
		} else {
			sid--;                      // delete vertex in stack
			p_idx--;

			last[qID[path[p_idx]]] = path[p_idx];
		}
	}

	//Memory cleanup
	Free(stack);
	Free(tree->p);
	Free(tree->dfsl);
	Free(tree->ptol);
	Free(tree);
	Free(last);
	Free(path);
	Free(tree_id);

	free_graph(g_dom);
	return z_array;


}


/*
 * Function: find_j
 * ----------------------------
 * It returns the largest j s.t. j < i and C[(j*ELEM_LIST_SIZE) +2 ] == 0
 *
 * @param C
 * @param i
 *
 * @returns
 *
 */
static VTYPE find_j(VTYPE *C, VTYPE i) {
	VTYPE k;
	k = (i * ELEM_LIST_SIZE); // I will start from the (i-1)th element and check the 3rd position (l_id)
	for ( k = (i * ELEM_LIST_SIZE); k != ELEM_LIST_SIZE; k -= ELEM_LIST_SIZE) {
		if (C[k - 1] == 0)
			return ( (k / ELEM_LIST_SIZE) - 1);
	}
	return 0;
}


/*
 * Function: list_element_compare
 * ---------------------------------------
 *  Auxiliary function, it compares two elements
 *
 *  @param A 		the first element
 *  @param B		the second element
 *
 *  @return		it returns -1 if A < B, 1 if A > B, 0 if A == B
 *
 */
static int list_element_compare(const void *A, const void *B) {
	VTYPE i;
	VTYPE *a, *b;

	a = (VTYPE*)A;
	b = (VTYPE*)B;
	for ( i = 0; i < ELEM_LIST_SIZE; i++) {
		if (a[i] < b[i]) {
			return -1;
		} else if ( a[i] > b[i] ) {
			return 1;
		}
	}
	return 0; // It should never reach this point
}


/*
 * Function: sort_list
 * -------------------------
 * qsort wrapper to sort a list
 *
 * @param list		the list to sort
 * @param n_elem	the number of elements in the list
 *
 * @return		the sorted list
 *
 */
static VTYPE *sort_list(VTYPE *list, VTYPE n_elem) {
	qsort(list, n_elem, ELEM_LIST_SIZE * sizeof(VTYPE), list_element_compare);
	return list;
}

/*
 * Function: find_w
 * ------------------------------------
 * For each vertices v, find w that is a child of d(h(v)) and ancestor of v in D and DR
 *
 * @param g		the input graph
 * @param d		the dominator array
 * @param h		the header array
 * @param dR		the dominator array of gR
 * @param dom_dfs	the DFS of the dominator tree of g
 * @param domR_dfs	the DFS of the dominator tree of gR
 * @param lastdescR	the last descendant array of gR
 *
 * @return
 *
 */
static VTYPE *find_w(GRAPH *g, VTYPE *d, VTYPE *h, VTYPE *dR, DFS_TREE *dom_dfs, DFS_TREE *domR_dfs, VTYPE *lastdescR) {

	VTYPE i, n, idxA, idxB;
	VTYPE *W;			// W= w s.t. child of d(h(v)) and ancestor of v in D and DR
	VTYPE *C; 	// Auxiliary lists
	VTYPE sizeA, sizeC;
	VTYPE j, v, w, pv, pw; // index of the largest index < i of C s.t. C[(j*3]+2] == 0, vertices v and w, preorder of v and w in DR
	n = g->n;

	// Allocate memory and initialize lists
	C = Malloc(2 * ELEM_LIST_SIZE * (n + 1) * sizeof(VTYPE));
	W = Malloc( (n + 1) *	 sizeof(VTYPE));

	memset(C, 0, 2 * ELEM_LIST_SIZE * (n + 1) * sizeof(VTYPE)); //The first ELEM_LIST_SIZE * n are used for B and the other for A
	memset(W, 0, (n + 1) * sizeof(VTYPE)); //The first ELEM_LIST_SIZE * n are used for B and the other for A

	// Construct list of pairs
	// and lists A and B
	idxB = ELEM_LIST_SIZE;
	idxA = (n + 1) * ELEM_LIST_SIZE;
	for ( i = 1; i <= n  ; i++, idxB += ELEM_LIST_SIZE) {
		VTYPE hv = h[i];
		if ( d[hv] != d[i] ) {
			C[idxA] 	= dom_dfs->dfsl[d[hv]]; // preorder number of (d(hv)) in D
			C[idxA + 1] = dom_dfs->dfsl[i];		// preorder number of i in D
			C[idxA + 2]	= 1;
			idxA += ELEM_LIST_SIZE;
		}
		C[idxB] 	= dom_dfs->dfsl[d[i]]; 		// preorder number of d[i] in D
		C[idxB + 1] = dom_dfs->dfsl[i];			// preorder number of i in D
		C[idxB + 2] = 0;
	}
	sizeA = (idxA / ELEM_LIST_SIZE);
	sizeC = sizeA;
	C = sort_list(C, sizeC);
	for ( i = (sizeC - 1) ; i != 0; i--) {
		VTYPE idxC = i * ELEM_LIST_SIZE;
		if ( 1 == C[idxC + 2]) {
			v = dom_dfs->ptol[ C[idxC + 1] ];	// v with preorder number C[idxC + 1] in D
			// Find the largest index j < i of list C s.t. C[(j*ELEM_LIST_SIZE) + 2] == 0
			j = find_j(C, i);
			w = dom_dfs->ptol[ C[(j * ELEM_LIST_SIZE) + 1] ]; 	// w with preorder number C[ (j*ELEM_LIST_SIZE)+1] ] in D

			pv = domR_dfs->dfsl[v];				// preorder number of v in DR
			pw = domR_dfs->dfsl[w];				// preorder number of w in DR
			if ( pw <= pv && pv <= lastdescR[pw]) {
				W[v] = w;
			}
		}
	}
	Free(C);
	return W;
}



/*
 * Function: compute_int_dfo
 * ----------------------------------------
 * PCD(v) = intersection of <D~(v), DR~(v)>
 *
 * @param g		the input graph
 * @param dom		the dominators
 * @param header	the headers
 * @param domR		the dominator of gR
 * @param headerR	the header of fT
 * @param dom_dfs	the DFS of the dominator tree of g
 * @param domR_dfs	the DFS of the dominator tree of gR
 * @param lastdescR	the last descendant array of gR
 * @param cscc		cscc
 *
 * @return 
 *
 */
static BUNDLES *compute_int_dfo(GRAPH *g, VTYPE *dom, VTYPE *header, VTYPE *domR,
		VTYPE *headerR, DFS_TREE *dom_dfs, DFS_TREE * domR_dfs, VTYPE *lastdescR, VTYPE *cscc) {

	// W <- Find_W(D,H,D',dsf_D, dfs_D')
	// Q <- Build_Q(D,D')
	// Z <- Find_Z(W,Q)

	BUNDLES *dfo;
	VTYPE i, n;
	VTYPE v, w, z,  c;
	VTYPE *W, *Q, *Z, *qID;
	LTYPE *sumC;
	LTYPE *sum_sqC;
	LTYPE *sizeC;
	LTYPE *startsumC, *endsumC, *startsum_sqC, *endsum_sqC, *startsize, *endsize;
	n = g->n;

	dfo 		 = Malloc( sizeof(BUNDLES));
	sumC 		 = Malloc( (n + 1) * sizeof(LTYPE));
	sum_sqC 	 = Malloc( (n + 1) * sizeof(LTYPE));
	startsumC 	 = Malloc( (n + 1) * sizeof(LTYPE));
	endsumC 	 = Malloc( (n + 1) * sizeof(LTYPE));
	startsum_sqC = Malloc( (n + 1) * sizeof(LTYPE));
	endsum_sqC   = Malloc( (n + 1) * sizeof(LTYPE));
	startsize 	 = Malloc( (n + 1) * sizeof(LTYPE));
	endsize 	 = Malloc( (n + 1) * sizeof(LTYPE));
	sizeC	 	 = Malloc( (n + 1) * sizeof(LTYPE));

	qID = Malloc( (n + 1) * sizeof(VTYPE));

	memset(sumC		, 0, (n + 1) * sizeof(LTYPE));
	memset(sum_sqC		, 0, (n + 1) * sizeof(LTYPE));
	memset(startsumC	, 0, (n + 1) * sizeof(LTYPE));
	memset(endsumC		, 0, (n + 1) * sizeof(LTYPE));
	memset(startsum_sqC	, 0, (n + 1) * sizeof(LTYPE));
	memset(endsum_sqC	, 0, (n + 1) * sizeof(LTYPE));
	memset(startsize	, 0, (n + 1) * sizeof(LTYPE));
	memset(endsize		, 0, (n + 1) * sizeof(LTYPE));
	memset(sizeC		, 0, (n + 1) * sizeof(LTYPE));

	memset(qID, 0, (n + 1) * sizeof(VTYPE));


	// Find_W
	W = find_w(g, dom, header, domR, dom_dfs, domR_dfs, lastdescR);

	// Build_Q
	Q = build_q(g, dom, domR, qID, dom_dfs);

	// Find_Z
	Z = find_z_dfs(g,W,qID, dom, g->root);
	// for each tuple <v,w,z> in >
	for ( i = 1 ; i < (n+1) ; i++) {
		LTYPE cscc_v = 0;
		LTYPE cscc_v_sq = 0;
		v = i;
		w = W[i];
		z = Z[i];
		cscc_v = cscc[v];
		cscc_v_sq = cscc_v * cscc_v;

		sizeC[z] += cscc_v;
		endsize[z] += cscc_v;
		startsize[w] 	+= cscc_v;
		if (cscc_v > 1) {
			startsumC[w] 	+= cscc_v;
			startsum_sqC[w] += cscc_v_sq;
			sumC[z] += cscc_v;
			endsumC[z] += cscc_v;
			sum_sqC[z] += cscc_v_sq;
			endsum_sqC[z] += cscc_v_sq;
		}
	}
	// Partial sum of CSCC for all common dominators i path z->w in Q
	// for each v in Q bottom up approach
	for ( i = n ; i != 0 ; i--) {
		c = dom_dfs->ptol[i];
		while ( (v = Q[c]) && (v != 0) ) {
			sumC[v] += sumC[c] - startsumC[c];
			sizeC[v] += sizeC[c] - startsize[c];
			sum_sqC[v] += sum_sqC[c] - startsum_sqC[c];
			Q[c] = 0; 		// set to 0 all the vertices already evaluated to avoid compute them multiple times
			c = v;

		}
	}

	Free(startsumC);
	Free(endsumC);
	Free(startsum_sqC);
	Free(endsum_sqC);
	Free(startsize);
	Free(endsize);
	Free(W);
	Free(Q);
	Free(qID);
	Free(Z);
	dfo->sumC	 = sumC;
	dfo->sum_sqC = sum_sqC;
	dfo->sizeC = sizeC;
	return dfo;
}
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
#define DEFAULT_ROOT 1
VTYPE set_root(VTYPE r1, VTYPE r2) {

	// If the root is not given from input (is still -1)
	if (r1 == -1) {
		// And the root is read from file (for example DIMACS)
		if (r2 != -1) {
			// Set the root read from file
			r1 = r2;
			// Set root to be equal to a default value
		} else {
			printf("Using default root value: %d\n", DEFAULT_ROOT);
			r1 = DEFAULT_ROOT;
		}
	}
	// Else we have nothing to do, just keep r1

	return r1;
}



/*
 * Function: get_graphs_from_cc
 * --------------------------------------------------
 * It computes the subgraphs of the connected components of a graph
 *
 * @param g 		the input graph
 * @param cc		the connected components
 *
 * @return		the subgraphs of cc of g
 *
 */
static GRAPH **get_graphs_from_cc(GRAPH *g, CC *cc) {

	/* Allocate and build the subgraphs of the Connected Components */
	GRAPH **sub_graphs;
	VTYPE v, w;						/* The vertex we are visting */
	VTYPE num_cc;
	VTYPE *vtos; 					/* The array containing the scc_id for each vertex */
	VTYPE *oton; 					/* The array containing the mapping between old and new index */
	VTYPE i, j, k;
	VTYPE *g_next_out, *g_next_in;
	num_cc = cc->number;



	sub_graphs = Calloc(num_cc, sizeof(GRAPH*));

	for (i = 0 ; i < num_cc; i++) {
		sub_graphs[i] = (GRAPH*) Calloc(1, sizeof(GRAPH));
	}

	for ( i = 0; i < num_cc; i++) {
		sub_graphs[i]->n = (cc->index[i + 1] - cc->index[i]);
		sub_graphs[i]->original_n = (g->original_n);
		sub_graphs[i]->g1out = Calloc( sub_graphs[i]->n + 2, sizeof(VTYPE));
		sub_graphs[i]->g1in = Calloc( sub_graphs[i]->n + 2, sizeof(VTYPE));
		sub_graphs[i]->labels = Malloc( (sub_graphs[i]->n + 1) * sizeof(VTYPE));
		if( g->removed_num > 0){
			sub_graphs[i]->removed = Malloc( (g->removed_num) * sizeof(VTYPE));
			memcpy(sub_graphs[i]->removed, g->removed, (g->removed_num) * sizeof(VTYPE));
		}
		sub_graphs[i]->removed_num = g->removed_num;
		sub_graphs[i]->root = 1; //static root
		sub_graphs[i]->cn_equiv = NULL; 
		sub_graphs[i]->cn_equiv = 0;
	}

	if ( 1 == num_cc) {
		// We will return a copy of g. It is not the best choice from performance point of view but it is safeier
		sub_graphs[0]->m = g->m;
		sub_graphs[0]->cn = g->cn;
		sub_graphs[0]->cn_dfo = g->cn_dfo;
		sub_graphs[0]->gin = Malloc(g->m * sizeof(VTYPE));
		sub_graphs[0]->gout = Malloc(g->m * sizeof(VTYPE));
		memcpy(sub_graphs[0]->gin, g->gin, g->m * sizeof(VTYPE));
		memcpy(sub_graphs[0]->gout, g->gout, g->m * sizeof(VTYPE));
		memcpy(sub_graphs[0]->g1in, g->g1in, (g->n + 2) * sizeof(VTYPE));
		memcpy(sub_graphs[0]->g1out, g->g1out, (g->n + 2) * sizeof(VTYPE));
		memcpy(sub_graphs[0]->labels, g->labels, (g->n + 1) * sizeof(VTYPE));
		sub_graphs[0]->root = g->root;
		return sub_graphs;
	}

	vtos = (VTYPE*) Malloc( (g->n + 1) * sizeof(VTYPE));
	oton = (VTYPE*) Malloc( (g->n + 1) * sizeof(VTYPE));
	for ( i = 0 ; i < num_cc; i++) {
		k = 1;
		for (j = cc->index[i]; j < cc->index[i + 1]; j++) { // For each vertex in CC <i>
			sub_graphs[i]->labels[k] = g->labels[cc->vertex[j]];
			vtos[cc->vertex[j]] = i + 1;
			oton[cc->vertex[j]] = k;
			k++;
		}
	}

	// Auxiliary array - the same for all the subgraphs
	g_next_out = Calloc((g->n + 2), sizeof(VTYPE));
	g_next_in = Calloc((g->n + 2), sizeof(VTYPE));

	for ( i = 0; i < num_cc; i++) {
		for (j = cc->index[i]; j < cc->index[i + 1]; j++) { /* For each vertex in CC <i> */
			v = cc->vertex[j];

			if ( 0 != g->g1out[v + 1]) {
				for ( k = g->g1out[v + 1] - 1; (k >= g->g1out[v]); k--) {
					/* check which neighbors of <v> are in the same component of <v> */
					w = g->gout[k];
					if ( vtos[w] == (i + 1) ) { /* <v> and <w> are in the same component */
						sub_graphs[i]->g1out[oton[v] + 1]++;
						sub_graphs[i]->g1in[oton[w] + 1]++;
					}
					if ( 0 == k) break; // k can be unsigned
				}
			}
		}
		for (j = 1; j <= sub_graphs[i]->n + 1 ; j++) {
			sub_graphs[i]->g1out[j]      += sub_graphs[i]->g1out[j - 1];
			g_next_out[j]  = sub_graphs[i]->g1out[j];
			sub_graphs[i]->g1in[j]       += sub_graphs[i]->g1in[j - 1];
			g_next_in[j]   = sub_graphs[i]->g1in[j];
		}
		sub_graphs[i]->m = sub_graphs[i]->g1out[sub_graphs[i]->n + 1];
		sub_graphs[i]->gout = Calloc( sub_graphs[i]->m +1, sizeof(VTYPE));
		sub_graphs[i]->gin = Calloc( sub_graphs[i]->m +1 , sizeof(VTYPE));


		for (j = cc->index[i]; j < cc->index[i + 1]; j++) {
			v = cc->vertex[j];
			if (g->g1out[v + 1] != 0) {
				for ( k = g->g1out[v + 1] - 1; k >= g->g1out[v] ; k--) {
					w = g->gout[k];
					if ( vtos[w] == (i + 1)) {
						VTYPE ov = oton[v];
						VTYPE ow = oton[w];
						sub_graphs[i]->gout[g_next_out[ov]] = ow;
						g_next_out[ov]++;
						sub_graphs[i]->gin[g_next_in[ow]] = ov;
						g_next_in[ow]++;
					}
					if ( 0 == k) break; // k can be unsigned
				}
			}
		}
	}


	Free(g_next_in);
	Free(g_next_out);
	Free(vtos);
	Free(oton);
	return sub_graphs;
}
/*
 * Function: free_bundles
 * ----------------------------
 *  It deallocates the memory belonging to a bundles
 *
 *  @param b		the input bundles data structure
 *
 *  @return
 *
 */
static void free_bundles(BUNDLES * b) {
	if ( NULL == b) {
		return;
	}

	if ( NULL != b->sumC) {
		Free(b->sumC);
	}

	if ( NULL != b->sum_sqC) {
		Free(b->sum_sqC);
	}
	if ( NULL != b->sizeC) {
		Free(b->sizeC);
	}

	Free(b);
}

/*
 * Function: free_dfs_tree
 * ----------------------------------
 *  It deallocates the memory belonging to a dfs tree
 *
 *  @param d		the input dfs tree
 *
 *  @return 
 *
 */
static void free_dfs_tree(DFS_TREE *d) {
	if ( NULL == d) {
		return;
	}

	if ( d->dfsl != NULL) {
		Free(d->dfsl);
	}

	if ( d->ptol != NULL) {
		Free(d->ptol);
	}

	if ( d->p != NULL) {
		Free(d->p);
	}

	Free(d);
}

/*
 * Function: free_clist
 * --------------------------------------
 *  It deallocates the memory belonging to a clist
 *
 *  @param c		the input clist
 *
 *  @return
 *
 */
static void free_clist(C_EDGE_LISTS *c) {
	if ( NULL == c) {
		return;
	}

	if ( NULL != c->cycle) {
		Free(c->cycle);
	}

	if ( NULL != c->last_cycle) {
		Free(c->last_cycle);
	}

	if ( NULL != c->cross) {
		Free(c->cross);
	}

	if ( NULL != c->last_cross) {
		Free(c->last_cross);
	}

	if ( NULL != c->next) {
		Free(c->next);
	}

	Free(c);
}


/*
 * Function: duplicate_graph
 * ------------------------------
 * It allows to duplicate a graph
 *
 * @param src		the input graph
 *
 * @return		the duplicated graph
 *
 */
static GRAPH * duplicate_graph(GRAPH * src){

	GRAPH *dst;
	VTYPE *gin, *g1in, *gout, *g1out, *labels, *removed, *cn_equiv;
	VTYPE n, m;

	n = src->n;
	m = src->m;
	// Allocate memory for the adjacency list
	dst   = (GRAPH*) Calloc(1	    , sizeof(GRAPH));
	g1in  = (VTYPE*) Malloc((n + 2) * sizeof(VTYPE));
	g1out = (VTYPE*) Malloc((n + 2) * sizeof(VTYPE));
	gin   = (VTYPE*) Malloc(    m   * sizeof(VTYPE));
	gout  = (VTYPE*) Malloc(    m   * sizeof(VTYPE));
	labels = (VTYPE*) Malloc((n + 1) * sizeof(VTYPE));
	if ( src->removed_num > 0 ){
		removed = (VTYPE*) Malloc((src->removed_num) * sizeof(VTYPE));
		memcpy(removed, src->removed, (src->removed_num) * sizeof(VTYPE));
	} else {
		removed = NULL;
	}
	if(src->cn_equiv_len > 0){
		cn_equiv = (VTYPE*) Malloc((src->cn_equiv_len) * sizeof(VTYPE));
	} else{
		cn_equiv = NULL;
	}

	memcpy(gin,     src->gin,                           m * sizeof(VTYPE));
	memcpy(gout,    src->gout,                          m * sizeof(VTYPE));
	memcpy(g1in,    src->g1in,                   (n + 2) * sizeof(VTYPE));
	memcpy(g1out,   src->g1out,                  (n + 2) * sizeof(VTYPE));
	memcpy(labels,  src->labels,                 (n + 1) * sizeof(VTYPE));
	memcpy(cn_equiv, src->cn_equiv,  (src->cn_equiv_len) * sizeof(VTYPE));

	dst->g1in = g1in;
	dst->g1out = g1out;
	dst->gin = gin;
	dst->gout = gout;

	dst->labels = labels;
	dst->n = n;
	dst->m = m;
	dst->root = src->root;
	dst->removed = removed;
	dst->removed_num = src->removed_num;
	dst->original_n = src->original_n; 
	dst->cn_equiv_len = src->cn_equiv_len;
	dst->cn_equiv = cn_equiv;
	dst->analyzed = src->analyzed;

	return dst;
}


/*
 * Function: add_subgraphs_to_pool
 * ----------------------------------
 *  It allows to add some sugraphs to an existing graph-pool
 *
 *  @param poool		the input graph pool
 *  @param subgraphs		the array of subgraphs to add
 *  @param num_subgraphs	the len of the array subgraphs
 *
 *  @return 			the input graph pool extended
 *
 */
static GRAPH_POOL * add_subgraphs_to_pool(GRAPH_POOL *pool, GRAPH **subgraphs, VTYPE num_subgraphs) {
	VTYPE i, j;
	VTYPE size;
	VTYPE best_pool;
	VTYPE idx_best_pool;

	if (NULL == pool) {
		FATAL_ERROR("pool is NULL");
	}

	size = pool->size;
	if ( pool->num_subgraphs > 0) {
		idx_best_pool = pool->idx_best;
		best_pool = pool->subgraphs[idx_best_pool]->cn_dfo;
	} else {
		best_pool = 0;
	}
	if ( (num_subgraphs + pool->num_subgraphs) >= size) {
		while ( size < (num_subgraphs + pool->num_subgraphs)) {
			size *= 2;
		}
		pool->subgraphs = Realloc(pool->subgraphs, (size * sizeof(GRAPH*)));
		pool->size = size;
	}
	//pool->start_index = pool->num_subgraphs;

	for ( i = 0, j = pool->num_subgraphs; i < num_subgraphs && j < size; i++, j++) {
		// removed for be complaint with bf output
		if (subgraphs[i]->n == 1) {
			j--;
			pool->isolated_node++;
			//    free_graph(subgraphs[i]);
			continue;
		}
		pool->subgraphs[j] = duplicate_graph(subgraphs[i]);
		if (subgraphs[i]->cn_dfo > best_pool) {
			best_pool = subgraphs[i]->cn_dfo;
			idx_best_pool = j;
		}
	}
	pool->num_subgraphs = j;

	return pool;
}

/*
 * Function: init_graph_pool
 * -------------------------------------------
 *  It initialises a graph_pool starting from a graph g
 *
 *  @param g			the input graph
 *  @param num_isolated_vertices
 *
 *  @return		the graph pool of g
 *
 */
static GRAPH_POOL * init_graph_pool(GRAPH *g, VTYPE *num_isolated_vertices) {
	GRAPH_POOL *pool;
	VTYPE i = 0;

	pool = Malloc(sizeof(GRAPH_POOL));

	pool->subgraphs = Calloc( INIT_POOL_SIZE , sizeof(GRAPH*));
	pool->size = INIT_POOL_SIZE;
	pool->num_subgraphs = 0;
	pool->start_index = pool->idx_best = 0;
	pool->idx_best_array = NULL;
	pool->idx_best_array_len = 0;

	INFO_MESSAGE("init_graph_pool - start computing scc on input graph");
	SCC *scc =  compute_scc_tarjan(g);


	INFO_MESSAGE("init_graph_pool - computing graphs from sccs");
	GRAPH **g_array_new = get_graphs_from_cc(g, scc);
#ifdef STATISTICS
	for ( i = 0 ; i < scc->number; i++) {
		fprintf(stderr, "[INFO]: The "TYPE_FORMAT"-th scc has "TYPE_FORMAT" vertices and "TYPE_FORMAT" edges\n", i, g_array_new[i]->n, g_array_new[i]->m);
	}
#endif
	*num_isolated_vertices = scc->number;
	pool->isolated_node = 0;
	add_subgraphs_to_pool(pool, g_array_new, scc->number);

	for(i = 0; i < scc->number; i++){
		free_graph(g_array_new[i]);
	}
	free_cc(scc);
	Free(g_array_new);
	return pool;
}


/*
 * Function: find_best_subgraph_in_pool
 * --------------------------------------
 *  It finds the subgraph that contains the best candidate vertex
 *
 *  @param pool		the input graph pool
 *
 * @return 		the index of the best subgraph
 *
 */ 
static VTYPE find_best_subgraph_in_pool(GRAPH_POOL *pool) {
	VTYPE i;
	VTYPE idx_best;
	LTYPE val_best;
	VTYPE scc_best;

	if (NULL == pool) {
		FATAL_ERROR("pool is NULL");
	}

	if ( pool->num_subgraphs < 1) {
		FATAL_ERROR("pool->num_subgraphs < 1");
	}

	if ( NULL == pool->subgraphs) {
		FATAL_ERROR("pool->subgraphs is NULL");
	}
	val_best = 0 ;
	idx_best = 0;
        scc_best = 0;
	for ( i = 0 ; i < pool->num_subgraphs; i++) {
		if (  (pool->subgraphs[i]->analyzed) && (pool->subgraphs[i]->cn_dfo > val_best)) {
		    idx_best = i;
		    scc_best = pool->subgraphs[i]->m;
		    val_best = pool->subgraphs[i]->cn_dfo;
		} else if (  (pool->subgraphs[i]->analyzed) && (pool->subgraphs[i]->cn_dfo == val_best)) {
		    if(pool->subgraphs[i]->m > scc_best){
		        idx_best = i;
		        scc_best = pool->subgraphs[i]->m;
		        val_best = pool->subgraphs[i]->cn_dfo;
                    }
                }
               
	}

	return idx_best;
}

/*
 * Function: get_and_remove_best_graph_from_pool
 * ---------------------------------------
 *  It looks for the subgraph in the pool that contains the best candidate to be removed.
 *  It removes this subgraph from the pool and it returns it
 *
 *  @param pool		the input graph pool
 *
 *  @return		the subgraph that contains the best candidate
 *
 */
static GRAPH * get_and_remove_best_graph_from_pool(GRAPH_POOL *pool) {
	if (NULL == pool) {
		FATAL_ERROR("pool is NULL");
	}

	if ( pool->num_subgraphs < 1) {
		FATAL_ERROR("pool->num_subgraphs < 1");
	}

	GRAPH * g = pool->subgraphs[pool->idx_best];

	if (pool->num_subgraphs > 1) {
		pool->subgraphs[pool->idx_best] = pool->subgraphs[pool->num_subgraphs - 1];
		pool->subgraphs[pool->num_subgraphs - 1 ] = NULL;
		pool->num_subgraphs -= 1;
		pool->idx_best = find_best_subgraph_in_pool(pool);
	} else {
		pool->idx_best = 0;
		pool->num_subgraphs -= 1;
	}
	return g;
}

/*
 * Function: free_graph_pool
 * -------------------------------
 *  It deallocates the memory belonging to a graph pool
 *
 *  @param p		the graph pool
 *
 *  @return
 *
 */
static void free_graph_pool(GRAPH_POOL *p) {
	VTYPE i;
	if ( NULL == p) {
		return;
	}
	for ( i = 0 ; i < p->num_subgraphs; i++) {
		if ( NULL != p->subgraphs[i]) {
			free_graph(p->subgraphs[i]);
		}
	}
	if ( NULL != p->subgraphs) {
		Free(p->subgraphs);
	}

	if(p->idx_best_array_len != 0 && p->idx_best_array != NULL){
		free(p->idx_best_array);
		p->idx_best_array_len = 0;
	}
	Free(p);
}


/*
 * Function: compute_cndp
 * -----------------------------------------
 * It computes the set of k critical nodes accordirg to CNH heuristic proposed in 
 * "Computing Critical Nodes in Directed Graphs. ACM Journal of Experimental Algorihmics 23 (2018)
 * N. Paudel, L. Georgiadis, G. F. Italiano"
 * 
 * @param g		the input graph in CSR format
 * @param k		the number of critical nodes to be removed
 * @param fout		the output file pointer (can be stdout)
 * @param ftime 	the output file pointer (can be stdout)
 *
 * @return		0 on success
 *
 */
  
int compute_cndp(GRAPH *g, VTYPE k, FILE* fout, FILE* ftime){

    GRAPH_POOL   *pool;
    GRAPH        *gR;                   // Reverse Graph
    DFS_TREE     *tree;                 // DFS tree (Input Graph)
    DFS_TREE     *treeR;                // DFS tree (Reverse Graph)
    C_EDGE_LISTS *clists;               // Auxiliary lists of cross and cycle edges (Input Graph)
    C_EDGE_LISTS *clistsR;              // Auxiliary lists of cross and cycle edges (Reverse Graph)
    VTYPE        *dom_tree;             // dominator tree (Input Graph)
    VTYPE        *dom_treeR;            // dominator tree (Reverse Graph)
    VTYPE        *ln_tree;              // loop nesting tree (Input Graph)
    VTYPE        *ln_treeR;             // loop nesting tree (Reverse Graph)
    VTYPE	 *cscc;			// cardinality of loop headed by each vertex (Input Graph)
    VTYPE	 *csccR;		// cardinality of loop headed by each vertex (Reverse Graph)

    
    VTYPE        num_isolated_node;     
    VTYPE        count_k;
    LTYPE        total_val;
    VTYPE        j;
    VTYPE        initial_scc;
    VTYPE	 n;
    VTYPE 	 m;
    VTYPE	 root;
    VTYPE	 root_degree;
    VTYPE 	 total_n;
    VTYPE 	 total_m;
    VTYPE	 num_rand_pick;
    

    // Others
    TIMER_DEF;

    // Initialize
    gR       = NULL;
    tree     = NULL;
    clists   = NULL;
    dom_tree = NULL;
    dom_treeR = NULL;
    ln_tree  = NULL;
    num_isolated_node = 0;
    num_rand_pick=0;
    
    n=g->n;
    m=g->m;
    root=g->root;
    
    if ( k > n) {
        fprintf(stdout, "k can not be greaten than n!\n");
        fprintf(stdout,"k will be n\n");
        k = n;
    }
    if (k == 0){
        k = n;
    }

    pool = init_graph_pool(g, &initial_scc);
    num_isolated_node =  initial_scc - pool->num_subgraphs;
    count_k = 1;
    j = 0;
    total_val = 0;
    if( initial_scc == 1){
        total_val = ((LTYPE)g->n * (g->n -1))/2;
    }
    else{

        for( j = 0, total_val = 0; j < pool->num_subgraphs; j++){
            LTYPE val = ((LTYPE)pool->subgraphs[j]->n * (pool->subgraphs[j]->n - 1))/2;
            total_val+= val;
        }
    }
    free_graph(g);
    fprintf(stdout,TXT_RED"STARTING CONNECTIVITY %llu\n"TXT_NO_COLOR,total_val);
    if(total_val == 1){
        //it is useless to compute the MCN
        fprintf(stdout,"It is useless compute MCN on G with connectivity 1\n");
        fprintf(stdout,"There are only two nodes connected, you may choose any of two as MCN\n");
        fprintf(stdout,"Exiting\n");
        return 0;
    }
    fprintf(fout,"%f,"TYPE_FORMAT","TYPE_FORMAT",%llu",0.0,(VTYPE)0,(VTYPE)0,total_val);
    fprintf(ftime,"Search&DeleteBest,Connectivity,OutputDump\n");
    total_n = 0;
    total_m = 0;
    for(j=0; j < pool->num_subgraphs; j++){
	    total_n += pool->subgraphs[j]->n;
	    total_m += pool->subgraphs[j]->m;
    }
    fprintf(fout,","TYPE_FORMAT","TYPE_FORMAT","TYPE_FORMAT,total_n,total_m, (pool->num_subgraphs + pool->isolated_node));

    fprintf(fout,","TYPE_FORMAT"\n",pool->isolated_node);
    TIMER_START;
    while (k > 0) {

        VTYPE i = 0;
        LTYPE best_dfo = 0;
        VTYPE idx_best = 0;
        for ( i = pool->start_index; i < pool->num_subgraphs; i++) {
            g = pool->subgraphs[i];
            if(g->analyzed != 0){
		if(g->cn_dfo > best_dfo){
                    best_dfo = g->cn_dfo;
                    idx_best = i;
		} else if(g->cn_dfo == best_dfo && g->m > pool->subgraphs[idx_best]->m){
                    idx_best = i;
		}
                continue;
            }
            n = g->n;
            m = g->m;

            root_degree = check_degree(root, g);
            if (root_degree < 0) {
                g->cn = g->cn_dfo = 0;
		g->analyzed=1;
                continue;
            }
            // Swap gin and gout to get the reverse graph from g
            gR = get_reverse_graph(g);

            // Compute the data structure needed to solve CNDP for g and gR
            // 0) Initialize data structures
            // 1) DFS + EDGE_LISTS_FOR_LOOP_NESTING (lists of cycle edges and cross edges)
            // 2) ISLT  (Simple Langauer Tarjan algorithm to find dominators, Iterative version)
            // 3) IHEAD (Tarjan-Georgiadis algorithm to find Headers of Loop Nesting Trees, Iterative version)

            // Step 0
            init_data_structures(&tree, &clists, n, m);
            if ((tree == NULL) || (clists == NULL)) {
                FATAL_ERROR("init_data_structures failed");
            }

            // First consider the graph g (then the reverse gR)
            // Step 1
            // Compute DFS with lists for Loop Nesting Tree
            // return DFS tree in tree and the lists in clists
            idfs_for_lnt(root, g, tree, clists);

            // Step 2
            // Compute dominators
            //
            // SHOULD BE POSSIBLE TO USE A DATA STRUCTURE ORDERED BY DFS AND NOT BY ORIGINAL VERTEX INDICES
            // I MEAN gin[k] and g1in[k] --> se k è l'ordine del dfs risparmio 2 accessi in memoria
            //
            dom_tree  = Malloc((n + 1) * sizeof(VTYPE));
            memset(dom_tree, 0, (n + 1)* sizeof(VTYPE));

            build_dom(g, tree, dom_tree);



            VTYPE *idom_tree = translate_dom(g, tree, root, dom_tree);

            // Step 3
            // Compute Loop Nesting tree
            // int build_header(GRAPH *g, DFS_TREE *tree, C_EDGE_LISTS *clists)

            ln_tree = build_header(g, tree, clists, &cscc);

            VTYPE *iln_tree = translate_header(g, tree, ln_tree);

            // Step 0 on gR


            init_data_structures(&treeR, &clistsR, n, m);
            if ((treeR == NULL) || (clistsR == NULL)) {
                FATAL_ERROR("init_data_structures failed");
            }


            // Step 1 on gR
            // Compute DFS with lists for Loop Nesting Tree
            // return DFS tree in tree and the lists in clists


            idfs_for_lnt(root, gR, treeR, clistsR);

            // Step 2 on gR

            dom_treeR  = Malloc((n + 1) * sizeof(VTYPE));
            memset(dom_treeR, 0, (n + 1)* sizeof(VTYPE));

            build_dom(gR, treeR, dom_treeR);

            VTYPE *idom_treeR = translate_dom(gR, treeR, root, dom_treeR);

            // Step 3 on gR
            ln_treeR = build_header(gR, treeR, clistsR, &csccR);

            VTYPE *iln_treeR = translate_header(gR, treeR, ln_treeR);
            VTYPE art_pts_num = 0;
            VTYPE *art_pts = NULL;
            art_pts = get_art_points(g, idom_tree, idom_treeR, &art_pts_num); 
            if( 0 == art_pts_num){
                WARNING("There are not articulation points\n");
                WARNING("Pick up the root as cndp. The removing of one vertex will create just one SCC\n");
                g->cn = g->root;
                g->cn_dfo = (g->n - 1) ; 
		g->randomly_chosen=1;
		if( g->cn_dfo > best_dfo){
		    best_dfo = g->cn_dfo;
		    idx_best = i;
		}
		g->analyzed=1;
                goto clean_no_art_pts;	
            }
            // Step 4 - Compute partial contributions to DFO
            // DFO-0  <- compute_DFO(D,H)
            // DFO-1  <- compute_DFO(D',H')
            // DFO-01 <- compute_INT_DFO(D,H,D',H')
            BUNDLES *dfo_0;
            GRAPH * g_dom = build_graph_from_array(g->n, g->n, idom_tree, root);
            VTYPE *lastdesc  = Malloc((n + 1) * sizeof(VTYPE));


            DFS_TREE *dfs_dom = simple_dfs(g->root, g_dom, lastdesc);

            VTYPE *icscc = translate_cscc(g, tree, cscc);


            dfo_0 = compute_dfo(g, idom_tree, iln_tree, icscc, dfs_dom->ptol);

            GRAPH * g_domR = build_graph_from_array(g->n, g->n, idom_treeR, root);
            VTYPE *lastdescR  = Malloc((n + 1) * sizeof(VTYPE));


            DFS_TREE *dfs_domR = simple_dfs(gR->root, g_domR, lastdescR);

            BUNDLES *dfo_1;
            VTYPE *icsccR = translate_cscc(gR, treeR, csccR);

            dfo_1 = compute_dfo(gR, idom_treeR, iln_treeR, icsccR, dfs_domR->ptol);

            BUNDLES * dfo_01 = compute_int_dfo(g, idom_tree, iln_tree, idom_treeR, iln_treeR, dfs_dom, dfs_domR, lastdescR, icscc);

            VTYPE w = 0;
            LTYPE out = 0;
            LTYPE dfo_tot = compute_sumDudNew(dfo_0, dfo_1, dfo_01, g, &w, &out,art_pts, art_pts_num);
            g->cn = w;
            g->cn_dfo = dfo_tot;
            if ( dfo_tot > best_dfo) {
                best_dfo = dfo_tot;
                idx_best = i;
            } else if (dfo_tot == best_dfo){
		if (g->m > pool->subgraphs[idx_best]->m){
		    idx_best=i;
		}	
	    }
            g->analyzed=1;
            free_bundles(dfo_0);
            free_bundles(dfo_1);
            free_bundles(dfo_01);

            free_graph(g_dom);
            free_graph(g_domR);

            Free(lastdesc);
            Free(lastdescR);

            Free(icscc);
            Free(icsccR);

            free_dfs_tree(dfs_dom);
            free_dfs_tree(dfs_domR);
clean_no_art_pts:
            free_clist(clists);
            free_clist(clistsR);

            free_dfs_tree(tree);
            free_dfs_tree(treeR);


            Free(dom_tree);
            Free(idom_tree);
            Free(dom_treeR);
            Free(idom_treeR);
            Free(cscc);
            Free(csccR);
            Free(ln_tree);
            Free(iln_tree);
            Free(ln_treeR);
            Free(iln_treeR);
            Free(gR); // Use free() on gR as it is only required for the struct
            Free(art_pts);
        }
        TIMER_START2;
        TIMER_START3;

        // We compare the old best with the new in order to identify the new critical node to remove
        if (pool-> num_subgraphs > 0 && pool->subgraphs[pool->idx_best]->cn_dfo < best_dfo ) {
            pool->idx_best = idx_best;
        } else if(pool-> num_subgraphs > 0 && pool->subgraphs[pool->idx_best]->cn_dfo == best_dfo){
	    if (pool->subgraphs[idx_best]->m > pool->subgraphs[pool->idx_best]->m){
		pool->idx_best=idx_best;
	    }
        } 


        g = get_and_remove_best_graph_from_pool(pool);
        total_val = total_val - g->cn_dfo; 
        TIMER_STOP3;
        float time_conn=(TIMER_ELAPSED3/1000000);
        // We evaluate all the subgraphs, now we need to choose which vertex will be removed and from which subgraph

        if(g->randomly_chosen){
	    num_rand_pick++;
	}

	TIMER_STOP;
	TIMER_START3;
    	fprintf(fout,"%f,"TYPE_FORMAT","TYPE_FORMAT",%llu",(TIMER_ELAPSED /  1000000.0 ),g->labels[g->cn],g->randomly_chosen,total_val);
        GRAPH *gn = remove_vertex_from_graph(g, g->cn);
      

        free_graph(g);

       
        SCC *scc =  compute_scc_tarjan(gn);
        
     

        GRAPH **g_array_new = get_graphs_from_cc(gn, scc);
        VTYPE j = 0;

        j = pool->num_subgraphs;
        add_subgraphs_to_pool(pool, g_array_new, scc->number);
	
	TIMER_STOP2;
        num_isolated_node += (scc->number - ( pool->num_subgraphs - j));

	total_n = 0;
	total_m = 0;
	for(j=0; j < pool->num_subgraphs; j++){
	    total_n += pool->subgraphs[j]->n;
	    total_m += pool->subgraphs[j]->m;
	}
	fprintf(fout,","TYPE_FORMAT","TYPE_FORMAT","TYPE_FORMAT,total_n,total_m, (pool->num_subgraphs + pool->isolated_node));
	
	fprintf(fout,","TYPE_FORMAT"\n",pool->isolated_node);
        TIMER_STOP3;
        fprintf(ftime,"%f,%f,%f\n",(TIMER_ELAPSED2/1000000.0),time_conn,(TIMER_ELAPSED3/1000000.0));

        if ( pool->num_subgraphs < 1) {
            WARNING("There are not scc with more than 1 vertex");
            WARNING("It is useless to continue to look for Critical Nodes- Every remaining vertices is a critical node - Stopping");
            k = 0;
        }

        VTYPE x = 0;
        for(x = 0; x < scc->number; x++){
        	free_graph(g_array_new[x]);
        }
        free(g_array_new);
        free_cc(scc);
        free_graph(gn);

        k--;
        count_k++;
	if(total_val == 0){
		WARNING("Connectivity is 0 - Stopping");
		k=0;
	}
    }
    fprintf(stdout, "[INFO] - Computation time = %f \n", TIMER_ELAPSED2 / 1000000.0);
    fprintf(stdout,"[INFO] - Randomly picked nodes (no articulation points) "TYPE_FORMAT"\n",num_rand_pick);
    free_graph_pool(pool);
    return 0;
}

