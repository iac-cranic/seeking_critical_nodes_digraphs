//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

/* File: utils.c */

#include "define.h"

void *Calloc(size_t nitems, size_t size) {

        void *ptr;

        ptr = calloc(nitems, size);
        if (!ptr) {
                fprintf(stderr, "Cannot allocate %zu items of size %zu\n", nitems, size);
                exit(EXIT_FAILURE);
        }
        return ptr;
}

void *Malloc(size_t sz) {
        void *ptr;

        ptr = malloc(sz);
        if (!ptr) {
                fprintf(stderr, "Cannot allocate %zu bytes\n", sz);
                exit(EXIT_FAILURE);
        }
        return ptr;
}

void *Realloc(void *ptr, size_t sz){
    void *lp;

    lp = (void *)realloc(ptr, sz);
    if (!lp && sz) {
        fprintf(stderr, "Cannot reallocate to %zu bytes...\n", sz);
        exit(EXIT_FAILURE);
    }
    return lp;
}

FILE *Fopen(const char *path, const char *mode)
{
    FILE *fp = NULL;
    fp = fopen(path, mode);
    if (!fp) {
        fprintf(stderr, "Cannot open file %s...\n", path);
        exit(EXIT_FAILURE);
    }
    return fp;
}

void Free(void *ptr){
	if(ptr){
		free(ptr);
	}
}

void usage(char *cmd)
{
    printf("usage: %s -f <input_file> -t <file_type> -o outfile [-r <root> -k <number_of_vertices_to_remove> (default =1) -h]\n", cmd);
    printf("supported file types:\n");
    printf("type 0 (DIMACS)\n");
    printf("type 1 (RMAT)\n");
    printf("\n-r specify the source vertex\n");
    printf("\n-k specify the number of vertex to remove. It should be 0 < k < #Vertices \n");
    printf("-h prints this short help\n");
	#ifdef LARGE
		printf("LARGE");
	#else
		printf("NOTDEFINED");
	#endif
}


char *Strdup(char *s){
	if( NULL == s){
		return s;
	}

	VTYPE len = strlen(s);
	char *d = Calloc(len+1, sizeof(char));
	memcpy(d,s,len);
	return d;
}


FILE * timeFileNameGen(char *outfilename){

   char  tfilen[MAX_BUFFER];
   if(-1 == snprintf(tfilen,MAX_BUFFER,"%s.time",outfilename)){
	fprintf(stderr,"[WARN] truncated output\n");
        exit(EXIT_FAILURE);
   }
   FILE * f=Fopen(tfilen,"w+");

   return f; 
}
