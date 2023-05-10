//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

#ifndef __CRANIC_BF_UTILITY_H_
#define __CRANIC_BF_UTILITY_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#ifndef DARWIN
  #include <bsd/string.h>
#endif

/* MACRO used to change the color of the standard output */
#define TXT_RED                     "\033[0;31m"
#define TXT_NO_COLOR                "\033[0m"
#define TXT_GREEN                   "\033[0;32m"
#define TXT_BLUE                    "\033[0;34m"
#define TXT_YELLOW                  "\033[0;33m"
#define TXT_CYAN                    "\033[0;36m"
#define TXT_PURPLE                  "\033[0;35m"

/* Some utility functions implemeted as MACRO */
#define MIN(a,b)  (a < b) ? a : b
#define MAX(a,b)  (a > b) ? a : b


	
/*
 * Function: Calloc
 * ------------------------------------------------
 * The function is a wrapper for the calloc function. The macro __FILE__ and 
 * __LINE__ are expected to be used as parameter for <file> and <line> 
 * respectivally.
 *
 * @param nitems	Number of items to allocate
 * @param size		Size of each item to be allocated
 * @param file     	File name from wich the function is called
 * @param line     	File's line from wich the funciton is called
 * @param msg     	Error message to show
 *
 * @return          A pointer to the memory allocated
 */
void *Calloc(size_t nitems, size_t size, const char *file, const unsigned long line, const char *msg);


/*
 * Function: Malloc
  ------------------------------------------------
 * The function is a wrapper for the malloc function. The macro __FILE__ and 
 * __LINE__ are expected to be used as parameters for <file> and <line> 
 * respectivally.
 *
 * @param size		Size of the memory to be allocated
 * @param file     	File name from wich the function is called
 * @param line     	File's line from wich the funciton is called
 * @param msg     	Error message to show
 *
 * @return          A pointer to the memory allocated
 */
void *Malloc(size_t size, const char *file, const unsigned long line, const char *msg);


/*
 * Function: Realloc
 * ------------------------------------------------
 * The function is a wrapper for the realloc function. The macro __FILE__ and 
 * __LINE__ are expected to be used as parameters for <file> and <line> 
 * respectivally.
 *
 * @param ptr		The pointer to the memory that has to be reallocated
 * @param size		The new size of the memory
 * @param file     	File name from wich the function is called
 * @param line     	File's line from wich the funciton is called
 * @param msg     	Error message to show
 *
 * @return          A pointer to the memory reallocated
 */
void *Realloc(void *ptr, size_t size, const char *file, const unsigned long line, const char *msg);


/*
 * Function: Fopen
 * ------------------------------------------------
 * The function is a wrapper for the fopen function. The macro __FILE__ and 
 * __LINE__ are expected to be used as parameters for <file> and <line> 
 * respectivally.
 *
 * @param path		The string containing the name of the file to be opened
 * @param mode		The string containing a file access mode
 * @param file     	File name from wich the function is called
 * @param line     	File's line from wich the funciton is called
 * @param msg     	Error message to show
 *
 * @return          A pointer to the file.
 */
FILE *Fopen(const char *path, const char *mode, const char *file, const unsigned long line, const char *msg);

#define TIME_SIZE 256
extern char time_str[TIME_SIZE];
/*
 * Function: get_time_string
 * ------------------------------------------------
 * The function returns the pointer to a string containing the current date in the
 * following format "%d-%m-%Y %H:%M:%S".
 *
 * @return          A pointer to a date string.
 */
char *get_time_string();


#define PATH_SIZE 512
extern char path[PATH_SIZE];
/*
 * Function: get_path
 * ------------------------------------------------
 * The function returns the pointer to a string containing the concatenation of a
 * dir path and a file name.
 *
 * @return          A pointer to a path string.
 */
const char *get_path(const char* dir, const char* file);

#endif
