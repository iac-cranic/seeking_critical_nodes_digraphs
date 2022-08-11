//This file belongs to the "Seeking Critical Nodes in Digraphs" project.
//The official repository is: "https://github.com/iac-cranic/seeking_critical_nodes_digraphs"
//This project is released under GPLv3; the full license file can be found in LICENSE file in the root of the repository.
//For any issue please contact us on github or at < cranic-info<AT>iac.rm.cnr.it >

/* File: timer.h */

#ifndef __CRANIC_CNH_TIMER_H_
#define __CRANIC_CNH_TIMER_H_

#include <sys/time.h>
#include <stdio.h>

#define TIMER_DEF     struct timeval temp_1, temp_2, temp_3, temp_4, temp_5, temp_6

#define TIMER_START   gettimeofday(&temp_1, (struct timezone*)0)
#define TIMER_START2   gettimeofday(&temp_3, (struct timezone*)0)
#define TIMER_START3   gettimeofday(&temp_5, (struct timezone*)0)

#define TIMER_STOP    gettimeofday(&temp_2, (struct timezone*)0)
#define TIMER_STOP2    gettimeofday(&temp_4, (struct timezone*)0)
#define TIMER_STOP3    gettimeofday(&temp_6, (struct timezone*)0)

#define TIMER_ELAPSED ((temp_2.tv_sec-temp_1.tv_sec)*1.e6+(temp_2.tv_usec-temp_1 .tv_usec))
#define TIMER_ELAPSED2 ((temp_4.tv_sec-temp_3.tv_sec)*1.e6+(temp_4.tv_usec-temp_3 .tv_usec))
#define TIMER_ELAPSED3 ((temp_6.tv_sec-temp_5.tv_sec)*1.e6+(temp_6.tv_usec-temp_5 .tv_usec))

#endif
