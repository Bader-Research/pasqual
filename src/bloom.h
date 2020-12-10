#ifndef __BLOOM_H__
#define __BLOOM_H__


#include "common.h"


// prototype of hash functions 
typedef INDEX_TYPE_T (*hashfunc_t) (char *, INDEX_TYPE_T);

// maximum number of hash functions
#define MAX_FUNC 5


int bloom_filter (unsigned char *array, INDEX_TYPE_T size, int nfunc,
                  char *str, char *r_str);


#endif /* __BLOOM_H__ */
