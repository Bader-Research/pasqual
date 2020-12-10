#ifndef __INDEX_H__
#define __INDEX_H__


#ifdef __OPENMP__
#include <omp.h>
#endif
#include "common.h"


// time of reading files
extern double time_importreads;


typedef struct _index_t
{
    // bloom filter
    unsigned char *bloom_array;
    INDEX_TYPE_T bloom_size;
    int n_hashfunc;

    // FM-Index
    INDEX_TYPE_T C[6];          // @ A C G T #
    INDEX_TYPE_T suffix_count[4];
    INDEX_TYPE_T prefix_count[4];
    u64 *OCC_bitmap[5];
    u64 *OCC_bitmap_bw[5];
    INDEX_TYPE_T *OCC_sample;
    INDEX_TYPE_T *OCC_sample_bw;

    // read index
    INDEX_TYPE_T *reads_index;
    INDEX_TYPE_T *reads_index_bw;

    // precomputed intervals
    INDEX_TYPE_T *interval_start;
    INDEX_TYPE_T *interval_end;
    u16 *interval_depth;
    INDEX_TYPE_T *interval_start_bw;
    INDEX_TYPE_T *interval_end_bw;
    u16 *interval_depth_bw;
    omp_lock_t *read_lock;

    // reads info
    char *reads;
    u16 read_len;
    INDEX_TYPE_T num_reads;
    u16 min_overlap;

    // SA buckets
    INDEX_TYPE_T size_bw[5];
    INDEX_TYPE_T offset_bw[5];
    INDEX_TYPE_T size_fw[5];
    INDEX_TYPE_T offset_fw[5];
    INDEX_TYPE_T max_size_fw;
    INDEX_TYPE_T max_size_bw;
} index_t;


/* -------------------------------------------------------------------------- */
/* definition of functions                                                    */

// build compressed index
index_t *construct_index (char *single_file,
                          char *pair_file1, char *pair_file2,
                          char *mate_file1, char *mate_file2,
                          int min_overlap, int read_len, int bloom_m,
                          int n_hashfunc);

// destroy index
void free_index (index_t * index);


#endif /* __INDEX_H__ */
