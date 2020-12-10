#include <stdio.h>
#include <math.h>
#ifdef __OPENMP__
#include <omp.h>
#endif
#ifdef __TIMING__
#include <time.h>
#include <sys/time.h>
#endif
#include <unistd.h>
#include <math.h>
#include <nmmintrin.h>

#include "sufsort.h"
#include "common.h"


/* -------------------------------------------------------------------------- */
/* global tables for triplet computation                                      */

char triplet_table[85] = { 0, 1, 22, 43, 64, 2, 7, 12, 17, 23,
    28, 33, 38, 44, 49, 54, 59, 65, 70, 75,
    80, 3, 4, 5, 6, 8, 9, 10, 11, 13,
    14, 15, 16, 18, 19, 20, 21, 24, 25, 26,
    27, 29, 30, 31, 32, 34, 35, 36, 37, 39,
    40, 41, 42, 45, 46, 47, 48, 50, 51, 52,
    53, 55, 56, 57, 58, 60, 61, 62, 63, 66,
    67, 68, 69, 71, 72, 73, 74, 76, 77, 78,
    79, 81, 82, 83, 84
};

char cool_table[SLOT_DEPTH][SLOT_NUM] =
    { {3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 0, 0, 0, 0, 0},
{18, 19, 20, 21, 24, 25, 26, 27, 29, 30, 31, 32, 0, 0, 0, 0, 0},
{34, 35, 36, 37, 39, 40, 41, 42, 45, 46, 47, 48, 50, 0, 0, 0, 0},
{51, 52, 53, 55, 56, 57, 58, 60, 61, 62, 63, 66, 67, 0, 0, 0, 0},
{68, 69, 71, 72, 73, 74, 76, 77, 78, 79, 81, 82, 83, 84, 0, 0, 0}
};

int cool[SLOT_DEPTH] = { 12, 12, 13, 13, 14 };


// SSE4 strcmp
inline int
STNI_strcmp (char *p1, int n1, char *p2, int n2)
{
#define MODE (_SIDD_SBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_LEAST_SIGNIFICANT  | _SIDD_NEGATIVE_POLARITY)

    __m128i smm1 = _mm_loadu_si128 ((__m128i *) p1);
    __m128i smm2 = _mm_loadu_si128 ((__m128i *) p2);
    int ResultIndex = 0;
    int remain_size = n1 < n2 ? n1 : n2;
    while (1)
    {
        ResultIndex = _mm_cmpistri (smm1, smm2, MODE);
        if (ResultIndex != 16 || ResultIndex > remain_size - 1)
        {
            break;
        }
        remain_size = remain_size - 16;
        p1 += 16;
        p2 += 16;
        smm1 = _mm_loadu_si128 ((__m128i *) p1);
        smm2 = _mm_loadu_si128 ((__m128i *) p2);
    }

    if (ResultIndex > remain_size - 1)
    {
        return n1 - n2;
    }
    p1 = (char *) &smm1;
    p2 = (char *) &smm2;

    return p1[ResultIndex] - p2[ResultIndex];

#undef MODE
}


void
transform_to_triplet_reverse (char *reads, int read_len,
                              INDEX_TYPE_T num_reads)
{
    INDEX_TYPE_T i;

#ifdef __OPENMP__
    #pragma omp parallel default(none) private(i) shared(triplet_table, read_len, reads, num_reads)
    {
#endif
        char *new_read = (char *) malloc (sizeof (char) * read_len);
#ifdef __OPENMP__
        #pragma omp for
#endif
        for (i = 0; i < num_reads; i++)
        {
            INDEX_TYPE_T start = read_len * (i + 1) - 1;
            int j;
            char t0 = reads[start] - 1;
            char t1 = reads[start - 1] - 1;
            char t2 = reads[start - 2] - 1;
            char s = t0 * 16 + t1 * 4 + t2;
            new_read[0] = reads[start] = triplet_table[s + 21];
            for (j = 1; j < read_len - 2; j++)
            {
                s = s & (char) 15;
                s = s * 4 + reads[start - j - 2] - 1;
                new_read[j] = reads[start - j] = triplet_table[s + 21];
            }
            s = s & (char) 15;
            new_read[read_len - 2] =
                reads[start - (read_len - 2)] = triplet_table[5 + s];
            s = s & (char) 3;
            new_read[read_len - 1] =
                reads[start - (read_len - 1)] = triplet_table[1 + s];
            memcpy (reads + i * read_len, new_read, sizeof (char) * read_len);
        }
        free (new_read);
#ifdef __OPENMP__
    }
#endif
}


void
transform_to_triplet (char *reads, int read_len, INDEX_TYPE_T num_reads)
{
    INDEX_TYPE_T i;
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < num_reads; i++)
    {
        INDEX_TYPE_T start = read_len * i;
        int j;
        char t0 = reads[start] - 1;
        char t1 = reads[start + 1] - 1;
        char t2 = reads[start + 2] - 1;
        char s = t0 * 16 + t1 * 4 + t2;
        reads[start] = triplet_table[s + 21];
        for (j = 1; j < read_len - 2; j++)
        {
            s = s & (char) 15;
            s = s * 4 + reads[start + j + 2] - 1;
            reads[start + j] = triplet_table[s + 21];
        }
        s = s & (char) 15;
        reads[start + read_len - 2] = triplet_table[5 + s];
        s = s & (char) 3;
        reads[start + read_len - 1] = triplet_table[1 + s];
    }
}


// when the bucket size is smaller than some threshold, use insertion sort
static void
insert_sort (INDEX_TYPE_T * sa, char *reads, INDEX_TYPE_T start, int pos,
             INDEX_TYPE_T num_elems, int read_len,
             INDEX_TYPE_T * interval_start, INDEX_TYPE_T * interval_end,
             u16 * interval_depth, omp_lock_t * read_lock,
             INDEX_TYPE_T bkt_offset, int min_overlap)
{
    int i;
    int j;
    INDEX_TYPE_T tmp;
    INDEX_TYPE_T *indices = &sa[start];
    INDEX_TYPE_T *new_index;
    INDEX_TYPE_T tmp_index;
    char *read_i;
    char *read_j;
    int len_i;
    int len_j;
    INDEX_TYPE_T read_num;

    new_index = (INDEX_TYPE_T *) malloc (sizeof (INDEX_TYPE_T) * num_elems);
    assert (new_index);

    tmp = indices[0];
    len_i = read_len - tmp % read_len - pos;
    if (len_i == 0 && pos <= min_overlap)
    {
        read_num = tmp / read_len;
        omp_set_lock (&read_lock[read_num]);
        if (interval_depth[read_num] < pos)
        {
            interval_start[read_num] = start + bkt_offset;
            interval_depth[read_num] = pos;
            interval_end[read_num] = num_elems;
        }
        omp_unset_lock (&read_lock[read_num]);
    }

    for (i = 1; i < num_elems; i++)
    {
        tmp = indices[i];
        tmp_index = i;
        read_i = &reads[tmp + pos];
        len_i = read_len - tmp % read_len - pos;
        if (len_i == 0 && pos <= min_overlap)
        {
            read_num = tmp / read_len;
            omp_set_lock (&read_lock[read_num]);
            if (interval_depth[read_num] < pos)
            {
                interval_start[read_num] = start + bkt_offset;
                interval_depth[read_num] = pos;
                interval_end[read_num] = num_elems;
            }
            omp_unset_lock (&read_lock[read_num]);
        }
        for (j = i - 1; j >= 0; j--)
        {
            len_j = read_len - indices[j] % read_len - pos;
            read_j = &reads[indices[j] + pos];
#ifdef __WITH_STNI__            
            if (STNI_strcmp (read_i, len_i, read_j, len_j) >= 0)
            {
                break;
            }
#else
            if (strncmp (read_i, read_j, min(len_i, len_j)) > 0)
            {
                break;
            }
            else if (strncmp (read_i, read_j, min(len_i, len_j)) == 0 && len_i >= len_j)
            {
                break;
            }
#endif
            indices[j + 1] = indices[j];
            new_index[j + 1] = new_index[j];
        }
        indices[j + 1] = tmp;
        new_index[j + 1] = tmp_index;
    }

    free (new_index);
}


// non-recursive sort
static void
process_slot_0 (INDEX_TYPE_T * X, char *reads, INDEX_TYPE_T start,
                INDEX_TYPE_T num_elems, int pos, int read_len,
                INDEX_TYPE_T * offsets, INDEX_TYPE_T * frequencies)
{
    INDEX_TYPE_T new_pos[NUM_SLOTS];
    int id;
    int read_of;
    INDEX_TYPE_T i;
    int bucket;
    INDEX_TYPE_T *indices = &X[start];
    INDEX_TYPE_T *new_indices =
        (INDEX_TYPE_T *) malloc (num_elems * sizeof (INDEX_TYPE_T));
    assert (new_indices != NULL);
    char *new_reads = (char *) malloc (num_elems * sizeof (char));
    assert (new_reads != NULL);
    memset (frequencies, 0, NUM_SLOTS * sizeof (INDEX_TYPE_T));

    for (i = 0; i < num_elems; ++i)
    {
        read_of = (indices[i] % read_len) + pos;
        id = (read_of >= read_len) ? (0) : (reads[indices[i] + pos]);
        new_reads[i] = id;
        frequencies[id]++;
    }

    new_pos[0] = offsets[0] = start;
    for (bucket = 1; bucket < NUM_SLOTS; ++bucket)
    {
        new_pos[bucket] = offsets[bucket] =
            offsets[bucket - 1] + frequencies[bucket - 1];
    }

    for (i = 0; i < num_elems; ++i)
    {
        id = new_reads[i];
        new_indices[new_pos[id] - start] = indices[i];
        new_pos[id]++;
    }

    memcpy (indices, new_indices, num_elems * sizeof (INDEX_TYPE_T));

    free (new_indices);
    free (new_reads);
}


// recursive sort
static void
process_slot_1 (INDEX_TYPE_T * X, INDEX_TYPE_T * new_indices, char *new_reads,
                char *reads, INDEX_TYPE_T start,
                INDEX_TYPE_T num_elems, int pos, int read_len,
                INDEX_TYPE_T * interval_start, INDEX_TYPE_T * interval_end,
                u16 * interval_depth, omp_lock_t * read_lock,
                INDEX_TYPE_T bkt_offset, int min_overlap)
{
    INDEX_TYPE_T offsets[NUM_SLOTS];
    INDEX_TYPE_T frequencies[NUM_SLOTS];
    int id;
    int read_of;
    INDEX_TYPE_T i;
    int bucket;
    INDEX_TYPE_T *indices = &X[start];
    int free_new_indices = 0;
    INDEX_TYPE_T read_num;

    if (new_indices == NULL)
    {
        new_indices =
            (INDEX_TYPE_T *) malloc (num_elems * sizeof (INDEX_TYPE_T));
        assert (new_indices != NULL);
        new_reads = (char *) malloc (num_elems * sizeof (char));
        assert (new_reads != NULL);
        free_new_indices = 1;
    }

    memset (frequencies, 0, NUM_SLOTS * sizeof (INDEX_TYPE_T));
    for (i = 0; i < num_elems; ++i)
    {
        read_of = (indices[i] % read_len) + pos;
        id = (read_of >= read_len) ? (0) : (reads[indices[i] + pos]);
        if (pos == read_len - (indices[i] % read_len) && pos <= min_overlap)
        {
            read_num = indices[i] / read_len;
            // precompute the interval for fast overlap finding
            omp_set_lock (&read_lock[read_num]);
            if (interval_depth[read_num] < pos)
            {
                interval_start[read_num] = start + bkt_offset;
                interval_depth[read_num] = pos;
                interval_end[read_num] = num_elems;
            }
            omp_unset_lock (&read_lock[read_num]);
        }
        new_reads[i] = id;
        frequencies[id]++;
    }

    offsets[0] = start;
    for (bucket = 1; bucket < NUM_SLOTS; ++bucket)
    {
        offsets[bucket] = offsets[bucket - 1] + frequencies[bucket - 1];
    }

    for (i = 0; i < num_elems; ++i)
    {
        id = new_reads[i];
        new_indices[offsets[id] - start] = indices[i];
        offsets[id]++;
    }

    memcpy (indices, new_indices, num_elems * sizeof (INDEX_TYPE_T));

    if (pos < read_len)
    {
        for (i = 21; i <= 84; i++)
        {
            int k = triplet_table[i];
            if (frequencies[k] > 1 &&
                frequencies[k] * (read_len - pos) < INSERT_SORT_SIZE)
            {
                insert_sort (X, reads, offsets[k - 1], pos + 3,
                             frequencies[k], read_len, interval_start,
                             interval_end, interval_depth, read_lock,
                             bkt_offset, min_overlap);
            }
            else if (frequencies[k] > 1)
            {
                process_slot_1 (X, new_indices, new_reads, reads,
                                offsets[k - 1], frequencies[k], pos + 3,
                                read_len, interval_start, interval_end,
                                interval_depth, read_lock, bkt_offset,
                                min_overlap);
            }
        }
    }

    if (free_new_indices == 1)
    {
        free (new_indices);
        free (new_reads);
    }
}


// NOTE: we do not malloc real suffix array in this step, only get the bucket size.
void
initial_sort (char *reads, INDEX_TYPE_T num_elems, bucket_t * bkt,
              INDEX_TYPE_T * bsize, INDEX_TYPE_T * boffset,
              INDEX_TYPE_T * max_size)
{
    INDEX_TYPE_T i;
    int j;
    int k;
    INDEX_TYPE_T *offsets = bkt->offs;
    INDEX_TYPE_T *frequencies = bkt->freqs;
    INDEX_TYPE_T msize;
    memset (frequencies, 0, NUM_SLOTS * sizeof (INDEX_TYPE_T));

#ifdef __OPENMP__
    omp_lock_t lock;
    omp_init_lock (&lock);
    #pragma omp parallel default(none) private(i, j, k) shared(lock, num_elems, reads, frequencies)
    {
#endif
        INDEX_TYPE_T freq[NUM_SLOTS];
        memset (freq, 0, NUM_SLOTS * sizeof (INDEX_TYPE_T));
#ifdef __OPENMP__
        #pragma omp for nowait
#endif
        for (i = 0; i < num_elems; i++)
        {
            freq[(int) reads[i]]++;
        }
#ifdef __OPENMP__
        omp_set_lock (&lock);
#endif
        for (j = 0; j < NUM_SLOTS; j++)
        {
            frequencies[j] += freq[j];
        }
#ifdef __OPENMP__
        omp_unset_lock (&lock);
        omp_destroy_lock (&lock);
    }
#endif

    msize = 0;
    for (k = 0; k < SLOT_DEPTH; k++)
    {
        offsets[k * SLOT_NUM] = 0;
        bsize[k] = frequencies[k * SLOT_NUM];
        for (j = 1; j < SLOT_NUM; j++)
        {
            bsize[k] += frequencies[k * SLOT_NUM + j];
            offsets[k * SLOT_NUM + j] =
                offsets[k * SLOT_NUM + j - 1] + frequencies[k * SLOT_NUM + j -
                                                            1];
        }
        if (bsize[k] > msize)
            msize = bsize[k];
    }

    boffset[0] = 0;
    for (k = 1; k < SLOT_DEPTH; k++)
    {
        boffset[k] = boffset[k - 1] + bsize[k - 1];
    }
    *max_size = msize;
}


void
suf_sort (INDEX_TYPE_T * sa, char *reads, int slot, INDEX_TYPE_T * new_idx,
          bucket_t * bkt0, INDEX_TYPE_T num_elems, int read_len,
          INDEX_TYPE_T * interval_start, INDEX_TYPE_T * interval_end,
          u16 * interval_depth, omp_lock_t * read_lock,
          INDEX_TYPE_T bkt_offset, int min_overlap)
{
    int size;
    int index;
    int new_index;
    bucket_t *bkt_queue = NULL;
    INDEX_TYPE_T i;
    int pos;
    int inter;
    int j;

    size = 0;
    for (j = 0; j < SEQ_DEPTH / 3; j++)
    {
        size += pow (64, j);
    }
    bkt_queue = (bucket_t *) calloc (size, sizeof (bucket_t));
    assert (bkt_queue != NULL);

    // STEP 1: postprocessing of intial sort
#ifdef __OPENMP__
    omp_lock_t *vlock =
        (omp_lock_t *) malloc (sizeof (omp_lock_t) * SLOT_NUM);
    assert (vlock != NULL);

    #pragma omp parallel default(none) private(i, j) shared(sa, slot, vlock, num_elems, reads, new_idx)
    {
        #pragma omp for
        for (j = 0; j < SLOT_NUM; j++)
        {
            omp_init_lock (&vlock[j]);
        }
#endif
        int jj;
        int k;
        INDEX_TYPE_T sa_buf[SLOT_NUM][128];
        INDEX_TYPE_T count[SLOT_NUM];
        memset (count, 0, SLOT_NUM * sizeof (INDEX_TYPE_T));
#ifdef __OPENMP__
        #pragma omp for
#endif
        for (i = 0; i < num_elems; i++)
        {
            INDEX_TYPE_T m;
            int n;
            int p;
            n = reads[i];
            if (n >= slot * SLOT_NUM && n < (slot + 1) * SLOT_NUM)
            {
                p = n - SLOT_NUM * slot;
                sa_buf[p][count[p]] = i;
                count[p]++;
                if (count[p] == 128)
                {
#ifdef __OPENMP__
                    omp_set_lock (&vlock[p]);
#endif
                    m = new_idx[n];
                    new_idx[n] = new_idx[n] + 128;
#ifdef __OPENMP__
                    omp_unset_lock (&vlock[p]);
#endif
                    memcpy (&sa[m], sa_buf[p], 128 * sizeof (INDEX_TYPE_T));
                    count[p] = 0;
                }
            }
        }
#ifdef __OPENMP__
        omp_set_lock (&vlock[0]);
#endif
        for (jj = 0; jj < SLOT_NUM; jj++)
        {
            for (k = 0; k < count[jj]; k++)
            {
                int p = slot * SLOT_NUM + jj;
                sa[new_idx[p]] = sa_buf[jj][k];
                new_idx[p] = new_idx[p] + 1;
            }
        }
#ifdef __OPENMP__
        omp_unset_lock (&vlock[0]);
    }
    #pragma omp parallel for
    for (j = 0; j < SLOT_NUM; j++)
    {
        omp_destroy_lock (&vlock[j]);
    }
    free (vlock);
#endif

    // STEP 2: parallel LS sort
    // TODO: Merge the code 
    inter = cool[slot];
#ifdef __OPENMP__
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (j = 0; j < inter; j++)
    {
        int k = cool_table[slot][j];
        bucket_t *curr = bkt0;
        if (curr->freqs[k] > 1)
        {
            bucket_t bkt;
            process_slot_0 (sa, reads, curr->offs[k],
                            curr->freqs[k], 3, read_len, bkt.offs, bkt.freqs);
            bkt_queue[j] = bkt;
        }
    }
    index = 0;
    new_index = inter;
    pos = 6;
    while (pos < SEQ_DEPTH)
    {
        inter = 64 * inter;
#ifdef __OPENMP__
        #pragma omp parallel for schedule(dynamic, 1)
#endif
        for (j = 0; j < inter; j++)
        {
            int m = j / 64;
            int n = j % 64;
            bucket_t *curr = &(bkt_queue[index + m]);
            int k = triplet_table[n + 21];
            if (curr->freqs[k] > 1)
            {
                bucket_t *bkt = &(bkt_queue[new_index + j]);
                process_slot_0 (sa, reads, curr->offs[k],
                                curr->freqs[k], pos,
                                read_len, bkt->offs, bkt->freqs);
            }
        }
        index = new_index;
        new_index += inter;
        pos += 3;
    }

    // STEP 3: recursive sort
#ifdef __OPENMP__
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (j = index * 64; j < new_index * 64; j++)
    {
        int m = j / 64;
        int n = j % 64 + 21;
        bucket_t *curr = &(bkt_queue[m]);
        int k = triplet_table[n];
        if (curr->freqs[k] > 1)
        {
            process_slot_1 (sa, NULL, NULL, reads, curr->offs[k],
                            curr->freqs[k], pos, read_len, interval_start,
                            interval_end, interval_depth, read_lock,
                            bkt_offset, min_overlap);
        }
    }

    free (bkt_queue);
}
