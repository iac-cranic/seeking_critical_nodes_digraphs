//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#include "utility.h"

char time_str[TIME_SIZE];
char path[PATH_SIZE];

void *Calloc(size_t nitems, size_t size, const char *file, const unsigned long line, const char *msg) {
	void *ptr;
    ptr = calloc(nitems, size);
    if (!ptr) {
		fprintf(stderr, "Function 'calloc' - Error msg: %s - FILE: %s LINE: %lu\n", msg, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}


void *Malloc(size_t size, const char *file, const unsigned long line, const char *msg) {
    void *ptr;
	ptr = malloc(size);  
	if (!ptr) {
        fprintf(stderr, "Function 'malloc' - Error msg: %s - FILE: %s LINE: %lu\n", msg, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}


void *Realloc(void *ptr, size_t size, const char *file, const unsigned long line, const char *msg){
    ptr = (void *)realloc(ptr, size);
    if (!ptr && size) {
        fprintf(stderr, "Function 'realloc' - Error msg: %s - FILE: %s LINE: %lu\n", msg, file, line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}


FILE *Fopen(const char *path, const char *mode, const char *file, const unsigned long line, const char *msg){
    FILE *fp = NULL;
    fp = fopen(path, mode);
    if (!fp) {
        fprintf(stderr, "Function 'fopen' - Error msg: %s - FILE: %s LINE: %lu\n", msg, file, line);
        exit(EXIT_FAILURE);
    }
    return fp;
}


char *get_time_string(){
	time_t rawtime;
	struct tm *info;
	
	time(&rawtime);
	info = localtime(&rawtime);
	
	strftime(time_str, sizeof(time_str), "%d-%m-%Y %H:%M:%S", info);
	
	return time_str;
}


const char *get_path(const char* dir, const char* file){
	if (strlcpy(path, dir, sizeof(path)) >= sizeof(path) ){
        fprintf(stderr, "(%s)[ERROR] The path %s. is too long, increase PATH_SIZE! Exit\n",get_time_string(), path);
		exit(EXIT_FAILURE);	
	}
	if (strlcat(path, "/", sizeof(path)) >= sizeof(path) ){
        fprintf(stderr, "(%s)[ERROR] The path %s. is too long, increase PATH_SIZE! Exit\n",get_time_string(), path);
		exit(EXIT_FAILURE);
	}
	if (strlcat(path, file, sizeof(path)) >= sizeof(path) ){
        fprintf(stderr, "(%s)[ERROR] The path %s. is too long, increase PATH_SIZE! Exit\n",get_time_string(), path);
		exit(EXIT_FAILURE);
	}
	return path;			  	
}

