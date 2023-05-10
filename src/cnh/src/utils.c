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
    printf("Usage: %s -f <INPUT_FILE> -t <FILE_TYPE> -o <OUTPUT_FILE> [-r <ROOT>] [-k <NUMBER>] [-h]\n", cmd);
    printf("\n\t-f <INPUT_FILE> \n");
    printf("\t   Read the graph from file <INPUT_FILE>.\n");
    printf("\t-t <FILE_TYPE> \n");
    printf("\t   The supported file types are 0 or 1; 0 for the DIMACS format and 1 for the RMAT format.\n");
    printf("\t-o <OUTPUT_FILE> \n");
    printf("\t   Write the results in the file <OUTPUT_FILE>.\n");
    printf("\t-r <ROOT> \n");
    printf("\t   Specify the source vertex\n");
    printf("\t-k <NUMBER>\n"); 
    printf("\t   Specify the number of vertex to remove. It should be 0 < NUMBER < #Vertices. The default value is 1. \n");
    printf("\t-h prints this short help\n\n");
//    #ifdef LARGE
//	printf("LARGE");
//    #else
//	printf("NOTDEFINED");
//    #endif
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
