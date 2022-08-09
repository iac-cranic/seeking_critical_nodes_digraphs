//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#ifndef __IGRAPH_CENTRALITIES_DEFINE_H_
#define __IGRAPH_CENTRALITIES_DEFINE_H_

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



#define MAX_BUFFER 4096

enum {
    DIMACS,            // 0
    RMAT,              // 1
} graph_type;
// OTHERS

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
 * Function: StringCopy
 * --------------------
 * Utility function to manage string copy with dinamically allocated memory
 *
 * @param s: the souce string to be copied
 *  *
 * @return: the pointer to dinamically allocated  memory containing a copy of s
 */
char *StringCopy(char *s);

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

#endif 
