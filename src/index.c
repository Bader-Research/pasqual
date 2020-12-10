#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "index.h"
#include "common.h"
#include "sufsort.h"
#include "bloom.h"
#ifdef __OPENMP__
#include <omp.h>
#endif
#ifdef __TIMING__
#include <time.h>
#include <sys/time.h>
#endif


double time_importreads = 0.0;

// reads buf size
#define READS_BUCK (1024 * 1024)


// build FM-Index
static void
build_OCC (index_t * index, INDEX_TYPE_T * sa, INDEX_TYPE_T * reads_index,
           INDEX_TYPE_T this_size, int *ridx, int dir,
           INDEX_TYPE_T this_start, int *last_bound, int *last_count)
{
    char *reads = index->reads;
    int mm = *ridx;
    int read_len = index->read_len;
    INDEX_TYPE_T i;
    int last = *last_bound;
    u64 **occ_bitmap;
    INDEX_TYPE_T occ_start;
    INDEX_TYPE_T *occ_sample;
    occ_start = this_start / 64;

    if (dir == 0)
    {
        occ_bitmap = index->OCC_bitmap;
        occ_sample = index->OCC_sample + this_start / 64 * 5;
    }
    else
    {
        occ_bitmap = index->OCC_bitmap_bw;
        occ_sample = index->OCC_sample_bw + this_start / 64 * 5;
    }

    // build BWT
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < this_size; i++)
    {
        int k = sa[i] % read_len;
        if (k != 0)
        {
            int nt = reads[sa[i] - 1];
            nt = (nt - 1) / 21 + 1;
            sa[i] = nt;
        }
        else
        {
            sa[i] = sa[i] / read_len + 20;
        }
    }
    for (i = 0; i < this_size; i++)
    {
        if (sa[i] >= 20)
        {
            reads_index[mm] = sa[i] - 20;
            mm++;
            sa[i] = 0;
        }
    }

    // build FM-Index
    for (i = 0; i < last; i++)
    {
        int nt;
        u64 m;
        nt = sa[i];
        last_count[nt]++;
        m = occ_bitmap[nt][i / 64 + occ_start];
        m = m | ((u64) 2 << (64 - last + i));
        occ_bitmap[nt][i / 64 + occ_start] = m;
    }
    occ_sample[0] = last_count[0];
    occ_sample[1] = last_count[1];
    occ_sample[2] = last_count[2];
    occ_sample[3] = last_count[3];
    occ_sample[4] = last_count[4];
    if (last != 0)
    {
        occ_start++;
        occ_sample += 5;
    }
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = last; i < this_size; i += 64)
    {
        int j;
        int nt;
        u64 m;
        int bound = (this_size - i) > 64 ? 64 : this_size - i;
        int st_count[5];
        st_count[0] = 0;
        st_count[1] = 0;
        st_count[2] = 0;
        st_count[3] = 0;
        st_count[4] = 0;
        for (j = 0; j < bound; j++)
        {
            nt = sa[i + j];
            st_count[nt]++;
            m = occ_bitmap[nt][i / 64 + occ_start];
            m = m | ((u64) 2 << j);
            occ_bitmap[nt][i / 64 + occ_start] = m;
        }
        if (this_size - i <= 64)
        {
            last_count[0] = st_count[0];
            last_count[1] = st_count[1];
            last_count[2] = st_count[2];
            last_count[3] = st_count[3];
            last_count[4] = st_count[4];
            *last_bound = 64 - bound;
        }
        occ_sample[0 + i / 64 * 5] = st_count[0];
        occ_sample[1 + i / 64 * 5] = st_count[1];
        occ_sample[2 + i / 64 * 5] = st_count[2];
        occ_sample[3 + i / 64 * 5] = st_count[3];
        occ_sample[4 + i / 64 * 5] = st_count[4];
    }

    *ridx = mm;
}


static int
import_singlereads (index_t * index, char *filename)
{
#ifdef __TIMING__
    struct timeval tv1, tv2;
    double time_pass;
    gettimeofday (&tv1, NULL);
#endif
    FILE *fp;

#ifdef DEBUG1
    FILE *fp_out;
    char tail[1000];
    char head[1000];
#endif

    char line[1000];
    int read_len;
    int i;
    INDEX_TYPE_T reads_count;
    INDEX_TYPE_T removed;
    INDEX_TYPE_T nt_count[4];
    char nt;
    char *r_line;
    char r_nt;
    INDEX_TYPE_T reads_size;
    int with_n;

    // init variables
    DPRINTF (3, "Importing reads from %s ... ", filename);
    read_len = (int) (index->read_len);
    nt_count[0] = 0;
    nt_count[1] = 0;
    nt_count[2] = 0;
    nt_count[3] = 0;
    if ((fp = fopen (filename, "r")) == NULL)
    {
        ERROR_EXIT ("open reads file %s failed.\n", filename);
    }

#ifdef DEBUG1
    if ((fp_out = fopen ("reads_filtered.fasta", "w+")) == NULL)
    {
        ERROR_EXIT ("open reads file %s failed.\n", filename);
    }
#endif

    r_line = (char *) malloc (sizeof (char) * (read_len + 1));
    assert (r_line != NULL);

    // count reads in file
    reads_count = 0;
    while (fgets (line, 1000, fp) != NULL)
    {
        if (line[0] == '>')
        {
            reads_count++;
        }
    }

    // import reads
    index->bloom_size = reads_count * index->bloom_size;
    index->bloom_array =
        calloc ((index->bloom_size + CHAR_BIT - 1) / CHAR_BIT,
                sizeof (unsigned char));
    assert (index->bloom_array != NULL);
    index->reads = (char *) malloc (sizeof (char) * read_len * READS_BUCK);
    assert (index->reads != NULL);
    reads_size = READS_BUCK;
    reads_count = 0;
    removed = 0;
    fseek (fp, 0, SEEK_SET);

    while (fgets (line, 1000, fp) != NULL)
    {
        if (line[0] != '>')
        {
#if defined(__INDEX_U32__)
            ERROR_EXIT ("incorrect reads file at line %d.\n",
                        reads_count * 2 + 2);
#elif defined(__INDEX_U64__)
            ERROR_EXIT ("incorrect reads file at line %lld.\n",
                        reads_count * 2 + 2);
#else
            ERROR_EXIT ("incorrect reads file at line %d.\n",
                        reads_count * 2 + 2);
#endif
        }
        else
        {
#ifdef DEBUG1
            memcpy (head, line, sizeof (char) * (strlen (line) - 1));
            head[strlen (line) - 1] = '\0';
#endif

            if (fgets (line, 1000, fp) == NULL)
            {
#if defined(__INDEX_U32__)
                ERROR_EXIT ("incorrect reads file at line %d.\n",
                            reads_count * 2 + 2);
#elif defined(__INDEX_U64__)
                ERROR_EXIT ("incorrect reads file at line %lld.\n",
                            reads_count * 2 + 2);
#else
                ERROR_EXIT ("incorrect reads file at line %d.\n",
                            reads_count * 2 + 2);
#endif
            }
            line[strlen (line) - 1] = '\0';

#ifdef DEBUG1
            memcpy (tail, line, sizeof (char) * strlen (line));
#endif

            if (strlen (line) != read_len)
            {
#if defined(__INDEX_U32__)
                ERROR_EXIT ("incorrect reads file at line %d.\n",
                            reads_count * 2 + 2);
#elif defined(__INDEX_U64__)
                ERROR_EXIT ("incorrect reads file at line %lld.\n",
                            reads_count * 2 + 2);
#else
                ERROR_EXIT ("incorrect reads file at line %d.\n",
                            reads_count * 2 + 2);
#endif
            }

            with_n = 0;
            // transform the symbols to integers
#ifdef __OPEN_MP__
        #pragma omp parallel for private(nt, r_nt)
#endif
            for (i = 0; i < read_len; i++)
            {
                nt = line[i];
                r_nt = 0;
                switch (nt)
                {
                case 'a':
                case 'A':
                    nt = 1;
                    r_nt = 4;
                    break;
                case 'c':
                case 'C':
                    nt = 2;
                    r_nt = 3;
                    break;
                case 'g':
                case 'G':
                    nt = 3;
                    r_nt = 2;
                    break;
                case 't':
                case 'T':
                    nt = 4;
                    r_nt = 1;
                    break;
                case 'n':
                case 'N':
                    with_n = 1;
                    break;
                default:
#if defined(__INDEX_U32__)
                    ERROR_EXIT ("incorrect reads file at line %d.\n",
                                reads_count * 2 + 2);
#elif defined(__INDEX_U64__)
                    ERROR_EXIT ("incorrect reads file at line %lld.\n",
                                reads_count * 2 + 2);
#else
                    ERROR_EXIT ("incorrect reads file at line %d.\n",
                                reads_count * 2 + 2);
#endif
                }
                line[i] = nt;
                r_line[read_len - 1 - i] = r_nt;
            }
            r_line[read_len] = '\0';

            // bloom filter duplicate reads
            if (with_n == 1 ||
                bloom_filter
                (index->bloom_array, index->bloom_size, index->n_hashfunc,
                 line, r_line) == 1)
            {
                removed++;
                continue;
            }

            // store unique reads
            nt = line[0];
            nt_count[nt - 1]++;
            index->prefix_count[nt - 1]++;
            nt = line[read_len - 1];
            nt_count[nt - 1]++;
            index->suffix_count[nt - 1]++;
            for (i = 1; i < read_len - 1; i++)
            {
                nt = line[i];
                nt_count[nt - 1]++;
            }
            memcpy (index->reads + reads_count * read_len, line,
                    sizeof (char) * read_len);
#ifdef DEBUG1
            fprintf (fp_out, ">READ %d %s\n", reads_count, head);
            fprintf (fp_out, "%s\n", tail);
#endif

            reads_count++;
            // extend the read buf
            if (reads_count == reads_size)
            {
                reads_size += READS_BUCK;
                index->reads =
                    (char *) realloc (index->reads,
                                      reads_size * read_len * sizeof (char));
                assert (index->reads != NULL);
            }
        }
    }
    index->num_reads = reads_count;
    index->C[0] = 0;
    index->C[1] = 0;
    index->C[2] = nt_count[0];
    index->C[3] = nt_count[0] + nt_count[1];
    index->C[4] = nt_count[0] + nt_count[1] + nt_count[2];
    index->C[5] = nt_count[0] + nt_count[1] + nt_count[2] + nt_count[3];

    DPRINTF (3, "done\n");
#if defined(__INDEX_U32__)
    DPRINTF (3,
             "  has imported %d reads (removed %d duplicate reads)\n",
             index->num_reads, removed);
#elif defined(__INDEX_U64__)
    DPRINTF (3,
             "  has imported %lld reads (removed %lld duplicate reads)\n",
             index->num_reads, removed);
#else
    DPRINTF (3,
             "  has imported %d (removed %d duplicate reads)\n",
             index->num_reads, removed);
#endif
    fclose (fp);
    free (index->bloom_array);

#ifdef DEBUG1
    fclose (fp_out);
#endif

#ifdef __TIMING__
    gettimeofday (&tv2, NULL);
    time_pass =
        (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec -
                                              tv1.tv_usec) / 1000.0;
    DPRINTF (3, "  takes %.3lf ms\n", time_pass);
    time_importreads = time_pass;
#endif

    return 0;
}


static void
build_idx (index_t * index, int dir)
{
    INDEX_TYPE_T i;
    int j;
    int mm;
    bucket_t bkt0;
    char *reads = index->reads;
    int read_len = (int) (index->read_len);
    INDEX_TYPE_T num_reads = index->num_reads;
    INDEX_TYPE_T num_elems = num_reads * read_len;
    INDEX_TYPE_T *sa;
    INDEX_TYPE_T *reads_index = NULL;
    INDEX_TYPE_T *max_size;
    INDEX_TYPE_T *bkt_size;
    INDEX_TYPE_T *bkt_offset;
    INDEX_TYPE_T new_index[NUM_SLOTS];
    INDEX_TYPE_T occ_rowlen;
    INDEX_TYPE_T *fix_count;
    int last_bound;
    int last_count[5];
    INDEX_TYPE_T s_count[5];
    INDEX_TYPE_T count[5];
    INDEX_TYPE_T *occ_sample;
    INDEX_TYPE_T *interval_start;
    INDEX_TYPE_T *interval_end;
    INDEX_TYPE_T *C = index->C;
    u16 *interval_depth;
    int min_overlap;
    u64 *occ_bitmap[5];
    occ_rowlen = (num_elems + 64) / 64;
    min_overlap = index->min_overlap;
#ifdef __OPENMP__
    omp_lock_t *read_lock = index->read_lock =
        (omp_lock_t *) malloc (index->num_reads * sizeof (omp_lock_t));
    assert (read_lock != NULL);
    for (i = 0; i < index->num_reads; i++)
    {
        omp_init_lock (&read_lock[i]);
    }
#endif

    // malloc OCC, read index
    if (dir == 0)               /* forward */
    {
        reads_index = index->reads_index =
            (INDEX_TYPE_T *) malloc (index->num_reads *
                                     sizeof (INDEX_TYPE_T));
        interval_start = index->interval_start =
            (INDEX_TYPE_T *) calloc (index->num_reads, sizeof (INDEX_TYPE_T));
        interval_end = index->interval_end =
            (INDEX_TYPE_T *) calloc (index->num_reads, sizeof (INDEX_TYPE_T));
        interval_depth = index->interval_depth =
            (u16 *) calloc (index->num_reads, sizeof (u16));
        assert (interval_start);
        assert (interval_end);
        assert (interval_depth);
        assert (index->reads_index != NULL);
        max_size = &(index->max_size_fw);
        bkt_size = index->size_fw;
        bkt_offset = index->offset_fw;
        fix_count = index->suffix_count;
        occ_bitmap[0] = index->OCC_bitmap[0] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[1] = index->OCC_bitmap[1] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[2] = index->OCC_bitmap[2] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[3] = index->OCC_bitmap[3] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[4] = index->OCC_bitmap[4] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        assert (occ_bitmap[0] != NULL || occ_bitmap[1] != NULL
                || occ_bitmap[2] != NULL || occ_bitmap[3] != NULL
                || occ_bitmap[4] != NULL);
        memset (occ_bitmap[0], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[1], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[2], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[3], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[4], 0, occ_rowlen * sizeof (u64));
        occ_sample = index->OCC_sample =
            (INDEX_TYPE_T *) malloc ((num_elems + 64) / 64 * 5 *
                                     sizeof (INDEX_TYPE_T));
        assert (occ_sample != NULL);
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < num_reads; i++)
        {
            char *read = reads + i * read_len;
            char ch = read[read_len - 1];
            interval_start[i] = C[(int) ch];
            interval_end[i] = C[(int) ch + 1] - C[(int) ch];
            interval_depth[i] = 1;
        }
    }
    else                        /* backward */
    {
        reads_index = index->reads_index_bw =
            (INDEX_TYPE_T *) malloc (index->num_reads *
                                     sizeof (INDEX_TYPE_T));
        interval_start = index->interval_start_bw =
            (INDEX_TYPE_T *) calloc (index->num_reads, sizeof (INDEX_TYPE_T));
        interval_end = index->interval_end_bw =
            (INDEX_TYPE_T *) calloc (index->num_reads, sizeof (INDEX_TYPE_T));
        interval_depth = index->interval_depth_bw =
            (u16 *) calloc (index->num_reads, sizeof (u16));
        assert (interval_start);
        assert (interval_end);
        assert (interval_depth);
        assert (index->reads_index_bw != NULL);
        max_size = &(index->max_size_bw);
        bkt_size = index->size_bw;
        bkt_offset = index->offset_bw;
        fix_count = index->prefix_count;
        occ_bitmap[0] = index->OCC_bitmap_bw[0] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[1] = index->OCC_bitmap_bw[1] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[2] = index->OCC_bitmap_bw[2] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[3] = index->OCC_bitmap_bw[3] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        occ_bitmap[4] = index->OCC_bitmap_bw[4] =
            (u64 *) malloc (occ_rowlen * sizeof (u64));
        assert (occ_bitmap[0] != NULL || occ_bitmap[1] != NULL
                || occ_bitmap[2] != NULL || occ_bitmap[3] != NULL
                || occ_bitmap[4] != NULL);
        memset (occ_bitmap[0], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[1], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[2], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[3], 0, occ_rowlen * sizeof (u64));
        memset (occ_bitmap[4], 0, occ_rowlen * sizeof (u64));
        occ_sample = index->OCC_sample_bw =
            (INDEX_TYPE_T *) malloc ((num_elems + 64) / 64 * 5 *
                                     sizeof (INDEX_TYPE_T));
        assert (occ_sample != NULL);
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < num_reads; i++)
        {
            char *read = reads + i * read_len;
            char ch = read[0];
            interval_start[i] = C[(int) ch];
            interval_end[i] = C[(int) ch + 1] - C[(int) ch];
            interval_depth[i] = 1;
        }
    }

    // tranform the reads to triplets
    if (dir == 0)
    {
        transform_to_triplet (reads, read_len, num_reads);
    }
    else
    {
        transform_to_triplet_reverse (reads, read_len, num_reads);
    }

    // intial sort
    initial_sort (reads, num_elems, &bkt0, bkt_size, bkt_offset, max_size);

    sa = (INDEX_TYPE_T *) malloc (*max_size * sizeof (INDEX_TYPE_T));
    assert (sa != NULL);
    mm = 0;
    last_bound = 0;
    memcpy (new_index, bkt0.offs, NUM_SLOTS * sizeof (INDEX_TYPE_T));

    // sort suffix array and build index in buckets
    for (j = 0; j < SLOT_DEPTH; j++)
    {
        suf_sort (sa, reads, j, new_index, &bkt0, num_elems, read_len,
                  interval_start, interval_end, interval_depth, read_lock,
                  bkt_offset[j], min_overlap);
        build_OCC (index, sa, reads_index, bkt_size[j], &mm, dir,
                   bkt_offset[j], &last_bound, last_count);
    }

    // prefix sum -> occ_sample
    count[0] = 0;
    count[1] = fix_count[0];
    count[2] = fix_count[1];
    count[3] = fix_count[2];
    count[4] = fix_count[3];
    for (i = 0; i < num_elems / 64; i++)
    {
        s_count[0] = occ_sample[(i + 1) * 5 + 0];
        s_count[1] = occ_sample[(i + 1) * 5 + 1];
        s_count[2] = occ_sample[(i + 1) * 5 + 2];
        s_count[3] = occ_sample[(i + 1) * 5 + 3];
        s_count[4] = occ_sample[(i + 1) * 5 + 4];
        occ_sample[(i + 1) * 5 + 0] = count[0] + occ_sample[i * 5 + 0];
        occ_sample[(i + 1) * 5 + 1] = count[1] + occ_sample[i * 5 + 1];
        occ_sample[(i + 1) * 5 + 2] = count[2] + occ_sample[i * 5 + 2];
        occ_sample[(i + 1) * 5 + 3] = count[3] + occ_sample[i * 5 + 3];
        occ_sample[(i + 1) * 5 + 4] = count[4] + occ_sample[i * 5 + 4];
        count[0] = s_count[0];
        count[1] = s_count[1];
        count[2] = s_count[2];
        count[3] = s_count[3];
        count[4] = s_count[4];
    }
    occ_sample[0 * 5 + 0] = 0;
    occ_sample[0 * 5 + 1] = fix_count[0];
    occ_sample[0 * 5 + 2] = fix_count[1];
    occ_sample[0 * 5 + 3] = fix_count[2];
    occ_sample[0 * 5 + 4] = fix_count[3];
    free (sa);

    for (i = 0; i < index->num_reads; i++)
    {
        omp_destroy_lock (&read_lock[i]);
    }
    free (read_lock);
}


static void
build_index (index_t * index)
{
#ifdef __TIMING__
    struct timeval tv1, tv2;
    double time_pass;
#endif
    INDEX_TYPE_T i;
    int read_len;
    INDEX_TYPE_T num_reads;
    char *reads;
    read_len = index->read_len;
    num_reads = index->num_reads;
    reads = index->reads;

    // forward index
    DPRINTF (3, "Sorting forward suffixes ... ");
#ifdef __TIMING__
    gettimeofday (&tv1, NULL);
#endif
    build_idx (index, 0);
    DPRINTF (3, "done\n");    
#ifdef __TIMING__
    gettimeofday (&tv2, NULL);
    time_pass =
        (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec -
                                              tv1.tv_usec) / 1000.0;
    DPRINTF (3, "  takes %.3lf ms\n", time_pass);
#endif

    // restore reads from triplets
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < read_len * num_reads; i++)
    {
        char nt = reads[i];
        nt = (nt - 1) / 21 + 1;
        reads[i] = nt;
    }

    // backward index
    DPRINTF (3, "Sorting backward suffixes ... ");
#ifdef __TIMING__
    gettimeofday (&tv1, NULL);
#endif
    build_idx (index, 1);
    DPRINTF (3, "done\n");
#ifdef __TIMING__
    gettimeofday (&tv2, NULL);
    time_pass =
        (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec -
                                              tv1.tv_usec) / 1000.0;
    DPRINTF (3, "  takes %.3lf ms\n", time_pass);
#endif

    // restore reads from triplets
#ifdef __OPENMP__
    #pragma omp parallel default(none) private(i) shared(read_len, num_reads, reads)
    {
#endif
        char *new_read = (char *) malloc (sizeof (char) * read_len);
#ifdef __OPENMP__
        #pragma omp for
#endif
        for (i = 0; i < num_reads; i++)
        {
            int j;
            for (j = 0; j < read_len; j++)
            {
                char nt = reads[i * read_len + j];
                nt = (nt - 1) / 21 + 1;
                new_read[read_len - 1 - j] = nt;
            }
            memcpy (reads + i * read_len, new_read, sizeof (char) * read_len);
        }
        free (new_read);
#ifdef __OPENMP__
    }
#endif
}


index_t *
construct_index (char *single_file,
                 char *pair_file1, char *pair_file2,
                 char *mate_file1, char *mate_file2,
                 int min_overlap, int read_len, int bloom_m, int n_hashfunc)
{
    index_t *index;
    index = (index_t *) malloc (sizeof (index_t));
    assert (index);

    memset (index, 0, sizeof (index_t));
    index->min_overlap = (u16) min_overlap;
    index->read_len = (u16) read_len;
    index->bloom_size = bloom_m;
    index->n_hashfunc = n_hashfunc;

    // import reads from file    
    if (single_file != NULL)
    {
        import_singlereads (index, single_file);
    }
    else if (pair_file1 != NULL && pair_file2 != NULL)
    {
        ERROR_EXIT
            ("We do not support paired end reads.\nPlease input them as single reads.\n");
    }
    else if (mate_file1 != NULL && mate_file2 != NULL)
    {
        ERROR_EXIT
            ("We do not support mate pair reads.\nPlease input them as single reads.\n");
    }

    // build compressed index
    build_index (index);

    return index;
}
