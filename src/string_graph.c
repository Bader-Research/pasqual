#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <emmintrin.h>

#include "string_graph.h"
#include "index.h"


inline INDEX_TYPE_T
OCCx (u64 ** OCC_bitmap, INDEX_TYPE_T * OCC_sample, int nt, INDEX_TYPE_T pos)
{
    INDEX_TYPE_T occ2 = 0;
    u64 occ_bitmap;
    int off;
    occ2 = OCC_sample[nt + pos / 64 * 5];
    off = pos % 64;
    occ_bitmap = OCC_bitmap[nt][pos / 64];
    occ_bitmap = occ_bitmap << (63 - off);
    occ2 += __builtin_popcountll (occ_bitmap);
    return occ2;
}


inline INDEX_TYPE_T
OCCx_LT (u64 ** OCC_bitmap, INDEX_TYPE_T * OCC_sample, int nt,
         INDEX_TYPE_T pos)
{
    INDEX_TYPE_T occ2 = 0;
    u64 occ_bitmap = 0;
    int off;
    int k;
    int pos0;
    pos0 = pos / 64;
    off = pos % 64;
    for (k = 0; k < nt; k++)
    {
        occ2 += OCC_sample[k + pos0 * 5];
        occ_bitmap |= OCC_bitmap[k][pos0];
    }
    occ_bitmap = occ_bitmap << (63 - off);
    occ2 += __builtin_popcountll (occ_bitmap);
    return occ2;
}


#define OCCfw(a, i) (OCCx(OCC_bitmap, OCC_sample, a, i))
#define OCCbw(a, i) (OCCx(OCC_bitmap_bw, OCC_sample_bw, a, i))
#define OCCfw_LT(a, i) (OCCx_LT(OCC_bitmap, OCC_sample, a, i))
#define OCCbw_LT(a, i) (OCCx_LT(OCC_bitmap_bw, OCC_sample_bw, a, i))


/* forward */
inline void
search_1 (INDEX_TYPE_T * C, u64 ** OCC_bitmap, INDEX_TYPE_T * OCC_sample,
          char ch, INDEX_TYPE_T * l, INDEX_TYPE_T * u)
{
    INDEX_TYPE_T occ2_l = 0;
    u64 occ_bitmap_l;
    int off_l;
    INDEX_TYPE_T occ2_u = 0;
    u64 occ_bitmap_u;
    int off_u;
    INDEX_TYPE_T vl = *l;
    INDEX_TYPE_T vu = *u + 1;

    occ2_l = OCC_sample[(int) ch + vl / 64 * 5];
    off_l = vl % 64;
    occ_bitmap_l = OCC_bitmap[(int) ch][vl / 64];
    occ_bitmap_l = occ_bitmap_l << (63 - off_l);
    occ2_l += __builtin_popcountll (occ_bitmap_l);

    occ2_u = OCC_sample[(int) ch + vu / 64 * 5];
    off_u = vu % 64;
    occ_bitmap_u = OCC_bitmap[(int) ch][vu / 64];
    occ_bitmap_u = occ_bitmap_u << (63 - off_u);
    occ2_u += __builtin_popcountll (occ_bitmap_u);

    *l = C[(int) ch] + occ2_l;
    *u = C[(int) ch] + occ2_u - 1;
}


/* forward and backward search */
inline void
search_2 (INDEX_TYPE_T * C, u64 ** OCC_bitmap, INDEX_TYPE_T * OCC_sample,
          char ch, INDEX_TYPE_T * l, INDEX_TYPE_T * u, INDEX_TYPE_T * ll,
          INDEX_TYPE_T * uu)
{
    INDEX_TYPE_T occ2_l = 0;
    u64 occ_bitmap_l;
    int off_l;
    INDEX_TYPE_T occ2_u = 0;
    u64 occ_bitmap_u;
    int off_u;
    INDEX_TYPE_T vl = *l;
    INDEX_TYPE_T vu = *u + 1;
    INDEX_TYPE_T occ2 = 0;
    u64 occ_bitmap_lt0 = 0;
    u64 occ_bitmap_lt1 = 0;
    int k;;

    off_l = vl % 64;
    off_u = vu % 64;
    for (k = 0; k < (int) ch; k++)
    {
        occ2 += OCC_sample[k + vu / 64 * 5];
        occ_bitmap_lt0 |= OCC_bitmap[k][vu / 64];
        occ2 -= OCC_sample[k + vl / 64 * 5];
        occ_bitmap_lt1 |= OCC_bitmap[k][vl / 64];
    }
    occ_bitmap_lt0 = occ_bitmap_lt0 << (63 - off_u);
    occ2 += __builtin_popcountll (occ_bitmap_lt0);
    occ_bitmap_lt1 = occ_bitmap_lt1 << (63 - off_l);
    occ2 -= __builtin_popcountll (occ_bitmap_lt1);
    *ll = *ll + occ2;

    occ2_l = OCC_sample[(int) ch + vl / 64 * 5];
    occ_bitmap_l = OCC_bitmap[(int) ch][vl / 64];
    occ_bitmap_l = occ_bitmap_l << (63 - off_l);
    occ2_l += __builtin_popcountll (occ_bitmap_l);
    occ2_u = OCC_sample[(int) ch + vu / 64 * 5];
    occ_bitmap_u = OCC_bitmap[(int) ch][vu / 64];
    occ_bitmap_u = occ_bitmap_u << (63 - off_u);
    occ2_u += __builtin_popcountll (occ_bitmap_u);

    *l = C[(int) ch] + occ2_l;
    *u = C[(int) ch] + occ2_u - 1;
    *uu = *ll + occ2_u - occ2_l - 1;
}


inline INDEX_TYPE_T
search_3 (INDEX_TYPE_T * C, u64 ** OCC_bitmap, INDEX_TYPE_T * OCC_sample,
          INDEX_TYPE_T vl, INDEX_TYPE_T vu)
{
    INDEX_TYPE_T occ2_l = 0;
    u64 occ_bitmap_l;
    int off_l;
    INDEX_TYPE_T occ2_u = 0;
    u64 occ_bitmap_u;
    int off_u;
    vu = vu + 1;

    occ2_l = OCC_sample[vl / 64 * 5];
    off_l = vl % 64;
    occ_bitmap_l = OCC_bitmap[0][vl / 64];
    occ_bitmap_l = occ_bitmap_l << (63 - off_l);
    occ2_l += __builtin_popcountll (occ_bitmap_l);

    occ2_u = OCC_sample[vu / 64 * 5];
    off_u = vu % 64;
    occ_bitmap_u = OCC_bitmap[0][vu / 64];
    occ_bitmap_u = occ_bitmap_u << (63 - off_u);
    occ2_u += __builtin_popcountll (occ_bitmap_u);

    return (occ2_u - occ2_l - 1);
}


// find overlaps
static void
find_intervals (INDEX_TYPE_T n, index_t * index, thread_resource_t * ts,
                int *overlaps_B, int *overlaps_F)
{
    int i;
    int k;
    INDEX_TYPE_T l;
    INDEX_TYPE_T u;
    INDEX_TYPE_T ll = 0;
    INDEX_TYPE_T uu = 0;
    char ch;
    char *read;
    int olp_B;
    int olp_F;
    int size_fw = 0;
    int size_bw = 0;
    interval_t *I_fw = ts->I_fw;
    interval_t *I_bw = ts->I_bw;;
    interval_t *I_next = NULL;
    INDEX_TYPE_T tmp;
    u64 **OCC_bitmap = index->OCC_bitmap;
    u64 **OCC_bitmap_bw = index->OCC_bitmap_bw;
    INDEX_TYPE_T *OCC_sample = index->OCC_sample;
    INDEX_TYPE_T *OCC_sample_bw = index->OCC_sample_bw;
    INDEX_TYPE_T *C = index->C;
    int read_len = (int) (index->read_len);
    int min_overlap = (int) (index->min_overlap);
    INDEX_TYPE_T *interval_start = index->interval_start;
    INDEX_TYPE_T *interval_end = index->interval_end;
    u16 *interval_depth = index->interval_depth;
    INDEX_TYPE_T *interval_start_bw = index->interval_start_bw;
    INDEX_TYPE_T *interval_end_bw = index->interval_end_bw;
    u16 *interval_depth_bw = index->interval_depth_bw;

    olp_B = olp_F = 0;

    char *reads = index->reads;
    read = reads + n * read_len;

    // find F-B intervals
    l = interval_start[n];
    u = l + interval_end[n] - 1;
    i = read_len - interval_depth[n] - 1;
    size_fw = 0;
    while (l < u && i > (read_len - min_overlap - 1))
    {
        ch = read[i];
        search_1 (C, OCC_bitmap, OCC_sample, ch, &l, &u);
        i = i - 1;
    }

    if (l < u)
    {
        ch = read[i + 1];
        ll = C[(int) ch];
        uu = C[(int) ch + 1] - 1;
        for (k = i + 2; k < read_len; k++)
        {
            ch = read[k];
            search_1 (C, OCC_bitmap_bw, OCC_sample_bw, ch, &ll, &uu);
        }
    }

    while (l < u && i >= 0)
    {
        tmp = search_3 (C, OCC_bitmap, OCC_sample, l, u);
        if (tmp >= 0)
        {
            I_next = &(I_fw[size_fw++]);
            assert (size_fw < read_len * 2);
            I_next->ll = ll;
            I_next->uu = ll + tmp;
            I_next->len = read_len - i - 1;
            I_next->type = FB;
            I_next->remain_size = tmp + 1;
            olp_F += tmp + 1;
        }
        ch = read[i];
        search_2 (C, OCC_bitmap, OCC_sample, ch, &l, &u, &ll, &uu);
        i = i - 1;
    }

    // find F-F intervals
    ch = ntc_table[(int) (read[read_len - 1])];
    l = C[(int) ch];
    u = C[(int) ch + 1] - 1;
    i = read_len - 2;

    while (l <= u && i > (read_len - min_overlap - 1))
    {
        ch = ntc_table[(int) (read[i])];
        search_1 (C, OCC_bitmap_bw, OCC_sample_bw, ch, &l, &u);
        i = i - 1;
    }
    if (l <= u)
    {
        ch = ntc_table[(int) read[i + 1]];
        ll = C[(int) ch];
        uu = C[(int) ch + 1] - 1;
        for (k = i + 2; k < read_len; k++)
        {
            ch = ntc_table[(int) read[k]];
            search_1 (C, OCC_bitmap, OCC_sample, ch, &ll, &uu);
        }
    }

    while (l <= u && i >= 0)
    {
        tmp = search_3 (C, OCC_bitmap_bw, OCC_sample_bw, l, u);
        if (tmp >= 0)
        {
            I_next = &(I_fw[size_fw++]);
            assert (size_fw < read_len * 2);
            I_next->ll = ll;
            I_next->uu = ll + tmp;
            I_next->len = read_len - i - 1;
            I_next->type = FF;
            I_next->remain_size = tmp + 1;
            olp_F += tmp + 1;
        }
        ch = ntc_table[(int) (read[i])];
        search_2 (C, OCC_bitmap_bw, OCC_sample_bw, ch, &l, &u, &ll, &uu);
        i = i - 1;
    }

    // find B-F intervals
    l = interval_start_bw[n];
    u = l + interval_end_bw[n] - 1;
    i = interval_depth_bw[n];
    size_bw = 0;
    while (l < u && i < min_overlap)
    {
        ch = read[i];
        search_1 (C, OCC_bitmap_bw, OCC_sample_bw, ch, &l, &u);
        i = i + 1;
    }
    if (l < u)
    {
        ch = read[i - 1];
        ll = C[(int) ch];
        uu = C[(int) ch + 1] - 1;
        for (k = i - 2; k >= 0; k--)
        {
            ch = read[k];
            search_1 (C, OCC_bitmap, OCC_sample, ch, &ll, &uu);
        }
    }

    while (l < u && i < read_len)
    {
        tmp = search_3 (C, OCC_bitmap_bw, OCC_sample_bw, l, u);
        if (tmp >= 0)
        {
            I_next = &(I_bw[size_bw++]);
            assert (size_bw < read_len * 2);
            I_next->ll = ll;
            I_next->uu = ll + tmp;
            I_next->len = i;
            I_next->type = BF;
            I_next->remain_size = tmp + 1;
            olp_B += tmp + 1;
        }
        ch = read[i];
        search_2 (C, OCC_bitmap_bw, OCC_sample_bw, ch, &l, &u, &ll, &uu);
        i = i + 1;
    }

    // find B-B intervals
    ch = ntc_table[(int) (read[0])];
    l = C[(int) ch];
    u = C[(int) ch + 1] - 1;
    i = 1;

    while (l <= u && i < min_overlap)
    {
        ch = ntc_table[(int) read[i]];
        search_1 (C, OCC_bitmap, OCC_sample, ch, &l, &u);
        i = i + 1;
    }
    if (l <= u)
    {
        ch = ntc_table[(int) read[i - 1]];
        ll = C[(int) ch];
        uu = C[(int) ch + 1] - 1;
        for (k = i - 2; k >= 0; k--)
        {
            ch = ntc_table[(int) read[k]];
            search_1 (C, OCC_bitmap_bw, OCC_sample_bw, ch, &ll, &uu);
        }
    }

    while (l <= u && i < read_len)
    {
        if (i >= min_overlap)
        {
            tmp = search_3 (C, OCC_bitmap, OCC_sample, l, u);
            if (tmp >= 0)
            {
                I_next = &(I_bw[size_bw++]);
                assert (size_bw < read_len * 2);
                I_next->ll = ll;
                I_next->uu = ll + tmp;
                I_next->len = i;
                I_next->type = BB;
                I_next->remain_size = tmp + 1;
                olp_B += tmp + 1;
            }
        }
        ch = ntc_table[(int) (read[i])];
        search_2 (C, OCC_bitmap, OCC_sample, ch, &l, &u, &ll, &uu);
        i = i + 1;
    }

    ts->size_fw = size_fw;
    ts->size_bw = size_bw;
    ts->max_fw = 0;
    ts->max_fw = 0;
    
    *overlaps_B = olp_B;
    *overlaps_F = olp_F;
}


// remove forward transitive edges
static void
extract_fw_edges (string_graph_t * graph, INDEX_TYPE_T n, index_t * index,
                  thread_resource_t * ts)
{
    int i;
    INDEX_TYPE_T k;
    INDEX_TYPE_T ll;
    INDEX_TYPE_T uu;    
    INDEX_TYPE_T lll;
    INDEX_TYPE_T uuu;
    int flag;
    int a;
    u16 max_len;
    int b;
    int size_Ia;
    interval_t *Ia;
    INDEX_TYPE_T overlap_read;
    int overlap_len;
    char type;
    int read_len = index->read_len;
    INDEX_TYPE_T *C = index->C;
    interval_t *I_F = ts->I_F;
    interval_t *I;
    short size_I;
    u64 **OCC_bitmap = index->OCC_bitmap;
    u64 **OCC_bitmap_bw = index->OCC_bitmap_bw;
    INDEX_TYPE_T *OCC_sample = index->OCC_sample;
    INDEX_TYPE_T *reads_index = index->reads_index;
    INDEX_TYPE_T *OCC_sample_bw = index->OCC_sample_bw;
    INDEX_TYPE_T *reads_index_bw = index->reads_index_bw;
   
    I = ts->I_fw + ts->pos_fw;
    size_I = ts->size_fw;
    flag = 0;
    overlap_len = read_len;
    type = FB;
    if (ts->max_fw + ts->level == read_len)
    {
        for (i = 0; i < size_I; i++)
        {
            ll = I[i].ll;
            uu = I[i].uu;
            overlap_len = I[i].len;
            type = I[i].type;
            if (overlap_len + ts->level != read_len)
                continue;
            if (type == FF)
            {
                lll = C[0] + OCCfw (0, ll);
                uuu = C[0] + OCCfw (0, uu + 1) - 1;
                //     search_1 (C, OCC_bitmap, OCC_sample, 0, &ll, &uu);
            }
            else
            {
                lll = C[0] + OCCbw (0, ll);
                uuu = C[0] + OCCbw (0, uu + 1) - 1;
                //    search_1 (C, OCC_bitmap_bw, OCC_sample_bw, 0, &ll, &uu);
            }

            if (lll <= uuu)
            {
                for (k = lll; k < uuu + 1; k++)
                {
                    if (I[i].type == FF)
                    {
                        overlap_read = reads_index[k];
                    }
                    else
                    {
                        overlap_read = reads_index_bw[k];
                    }
                    flag = 1;
                    if (overlap_read != n)
                    {                       
                        I_F[ts->size_I_F].ll = overlap_read;
                        I_F[ts->size_I_F].remain_size = ts->max_trans;
                        I_F[ts->size_I_F].len = overlap_len;
                        I_F[ts->size_I_F].type = type;
                        if (I_F[1].uu < ts->max_trans)
                        {
                            I_F[1].uu = ts->max_trans;
                        }
                        if (ts->size_I_F + 1 > 65535)
                        {
                            ERROR_EXIT ("assembly failed. Try a larger min overlap length.\n");
                        }
                        ts->size_I_F = ts->size_I_F + 1;
                        if (ts->size_I_F == I_F[0].uu)
                        {
                            I_F[0].uu = I_F[0].uu + read_len;
                            ts->I_F = I_F =
                                (interval_t *) realloc (I_F,
                                                        sizeof (interval_t) *
                                                        I_F[0].uu);
                            assert (I_F != NULL);
                        }
#if defined(__INDEX_U32__)
                        DPRINTF (5,
                                 "read %d add overlap type %d to read %d with len %d\n",
                                 n, type, overlap_read, overlap_len);
#elif defined(__INDEX_U64__)
                        DPRINTF (5,
                                 "read %lld add overlap type %d to read %lld with len %d\n",
                                 n, type, overlap_read, overlap_len);
#else
                        DPRINTF (5,
                                 "read %d add overlap type %d to read %d with len %d\n",
                                 n, type, overlap_read, overlap_len);
#endif
                    }                   
                }
                break;
            }
        }

        if (flag == 1)
        {
            return;
        }
    }

    Ia = I + size_I;
    ts->pos_fw = ts->pos_fw + size_I;
    for (a = 1; a < 5; a++)
    {
        size_Ia = 0;
        b = ntc_table[a];
        max_len = 0;
        ts->max_trans = 0; 
        for (i = 0; i < size_I; i++)
        {
            if (I[i].remain_size > 0)
            {
                ll = I[i].ll;
                uu = I[i].uu;
                if (I[i].type == FF)
                {
                    search_1 (C, OCC_bitmap, OCC_sample, b, &ll, &uu);
                }
                else
                {
                    search_1 (C, OCC_bitmap_bw, OCC_sample_bw, a, &ll, &uu);
                }
                I[i].remain_size = I[i].remain_size - (uu - ll + 1);
                if (ll <= uu)
                {
                    Ia[size_Ia].ll = ll;
                    Ia[size_Ia].uu = uu;
                    Ia[size_Ia].len = I[i].len;
                    Ia[size_Ia].type = I[i].type;
                    Ia[size_Ia].remain_size = uu - ll + 1;
                    ts->max_trans += uu - ll + 1;
                    if (max_len < I[i].len)
                        max_len = I[i].len;
                    size_Ia++;
                    assert (size_Ia + ts->size_fw + ts->pos_fw <
                            read_len * 2 * read_len * 2);
                }
            }
        }
        if (size_Ia == 0)
            continue;
        ts->size_fw = size_Ia;
        ts->level++;
        ts->max_fw = max_len;
        extract_fw_edges (graph, n, index, ts);
        ts->level--;
    }
    ts->pos_fw = ts->pos_fw - size_I;
}


// remove backward transitive edges
static void
extract_bw_edges (string_graph_t * graph, INDEX_TYPE_T n, index_t * index,
                  thread_resource_t * ts)
{
    int i;
    INDEX_TYPE_T k;
    INDEX_TYPE_T ll;
    INDEX_TYPE_T uu;
    INDEX_TYPE_T lll;
    INDEX_TYPE_T uuu;    
    int flag;
    int a;
    u16 max_len;
    int b;
    int size_Ia;
    interval_t *Ia;
    INDEX_TYPE_T overlap_read;
    int overlap_len;
    char type;
    int read_len = index->read_len;
    INDEX_TYPE_T *C = index->C;
    interval_t *I_B = ts->I_B;
    interval_t *I = ts->I_bw + ts->pos_bw;
    short size_I = ts->size_bw;
    u64 **OCC_bitmap = index->OCC_bitmap;
    u64 **OCC_bitmap_bw = index->OCC_bitmap_bw;
    INDEX_TYPE_T *OCC_sample = index->OCC_sample;
    INDEX_TYPE_T *reads_index = index->reads_index;
    INDEX_TYPE_T *OCC_sample_bw = index->OCC_sample_bw;
    INDEX_TYPE_T *reads_index_bw = index->reads_index_bw;

    I = ts->I_bw + ts->pos_bw;
    size_I = ts->size_bw;
    flag = 0;
    overlap_len = read_len;
    type = BB;
    if (ts->max_bw + ts->level == read_len)
    {
        for (i = 0; i < size_I; i++)
        {
            ll = I[i].ll;
            uu = I[i].uu;
            overlap_len = I[i].len;
            type = I[i].type;
            if (overlap_len + ts->level != read_len)
                continue;
            if (I[i].type == BB)
            {
                lll = C[0] + OCCbw (0, ll);
                uuu = C[0] + OCCbw (0, uu + 1) - 1;
                //     search_1 (C, OCC_bitmap_bw, OCC_sample_bw, 0, &ll, &uu);
            }
            else
            {
                lll = C[0] + OCCfw (0, ll);
                uuu = C[0] + OCCfw (0, uu + 1) - 1;
                //     search_1 (C, OCC_bitmap, OCC_sample, 0, &ll, &uu);
            }
            if (lll <= uuu)
            {
                overlap_len = I[i].len;
                type = I[i].type;
                for (k = lll; k < uuu + 1; k++)
                {
                    if (I[i].type == BB)
                    {
                        overlap_read = reads_index_bw[k];
                    }
                    else
                    {
                        overlap_read = reads_index[k];
                    }
                    flag = 1;
                    if (overlap_read != n)
                    {                       
                        I_B[ts->size_I_B].ll = overlap_read;
                        I_B[ts->size_I_B].remain_size = ts->max_trans;
                        I_B[ts->size_I_B].len = overlap_len;
                        I_B[ts->size_I_B].type = type;
                        if (I_B[1].uu < ts->max_trans)
                        {
                            I_B[1].uu = ts->max_trans;
                        }
                        if (ts->size_I_B + 1 > 65535)
                        {
                            ERROR_EXIT ("assembly failed. Try a larger min overlap length.\n");
                        } 
                        ts->size_I_B = ts->size_I_B + 1;                                  
                        if (ts->size_I_B == I_B[0].uu)
                        {
                            I_B[0].uu = I_B[0].uu + read_len;
                            ts->I_B = I_B =
                                (interval_t *) realloc (I_B,
                                                        sizeof (interval_t) *
                                                        I_B[0].uu);
                            assert (I_B != NULL);
                        }                     
#if defined(__INDEX_U32__)
                        DPRINTF (5,
                                 "read %d add overlap type %d to read %d with len %d\n",
                                 n, type, overlap_read, overlap_len);
#elif defined(__INDEX_U64__)
                        DPRINTF (5,
                                 "read %lld add overlap type %d to read %lld with len %d\n",
                                 n, type, overlap_read, overlap_len);
#else
                        DPRINTF (5,
                                 "read %d add overlap type %d to read %d with len %d\n",
                                 n, type, overlap_read, overlap_len);
#endif
                    }                   
                }
                break;
            }
        }

        if (flag == 1)
        {
            return;
        }
    }

    Ia = I + size_I;
    ts->pos_bw = ts->pos_bw + size_I;
    for (a = 1; a < 5; a++)
    {
        size_Ia = 0;
        b = ntc_table[a];
        max_len = 0;
        ts->max_trans = 0;
        for (i = 0; i < size_I; i++)
        {
            if (I[i].remain_size > 0)
            {
                ll = I[i].ll;
                uu = I[i].uu;
                if (I[i].type == BB)
                {
                    search_1 (C, OCC_bitmap_bw, OCC_sample_bw, a, &ll, &uu);
                }
                else
                {
                    search_1 (C, OCC_bitmap, OCC_sample, b, &ll, &uu);
                }
                I[i].remain_size = I[i].remain_size - (uu - ll + 1);

                if (ll <= uu)
                {
                    Ia[size_Ia].ll = ll;
                    Ia[size_Ia].uu = uu;
                    Ia[size_Ia].len = I[i].len;
                    Ia[size_Ia].type = I[i].type;
                    Ia[size_Ia].remain_size = uu - ll + 1;
                    ts->max_trans += uu - ll + 1;
                    if (max_len < I[i].len)
                        max_len = I[i].len;
                    size_Ia++;
                    assert (size_Ia + ts->size_fw + ts->pos_fw <
                            read_len * 2 * read_len * 2);
                }
            }
        }
        if (size_Ia == 0)
            continue;
        ts->size_bw = size_Ia;
        ts->level++;
        ts->max_bw = max_len;
        extract_bw_edges (graph, n, index, ts);
        ts->level--;
    }
    ts->pos_bw = ts->pos_bw - size_I;
}


string_graph_t *
construct_string_graph (index_t * index)
{
#ifdef __TIMING__
    struct timeval tv1, tv2;
    double time_pass;
    gettimeofday (&tv1, NULL);
#endif
    DPRINTF (3, "Building overlap graph ... ");
    string_graph_t *graph;
    double overlap_len;
    INDEX_TYPE_T overlaps;
    graph = (string_graph_t *) malloc (sizeof (string_graph_t));
    assert (graph != NULL);
    graph->nodes = (node_t *) malloc (index->num_reads * sizeof (node_t));
    assert (graph->nodes);
    memset (graph->nodes, 0, index->num_reads * sizeof (node_t));
    graph->num_reads = index->num_reads;
    graph->read_len = index->read_len;
    graph->min_overlap = index->min_overlap;
    overlaps = 0;
    overlap_len = 0;

#ifdef DEBUG1
    FILE *fp;
    FILE *fp1;
    fp = fopen ("overlaps.dat", "w+");
    fp1 = fopen ("len.dat", "w+");
#endif

#ifndef DEBUG1
#ifdef __OPENMP__
    #pragma omp parallel default(none) shared(index, graph, stderr, overlaps, overlap_len)
    {
#endif
#endif
        thread_resource_t *ts =
            (thread_resource_t *) malloc (sizeof (thread_resource_t));
        assert (ts != NULL);
        ts->I_fw =
            (interval_t *) malloc (sizeof (interval_t) * index->read_len * 2 *
                                   index->read_len * 2);
        ts->I_bw =
            (interval_t *) malloc (sizeof (interval_t) * index->read_len * 2 *
                                   index->read_len * 2);
        ts->I_F =
            (interval_t *) malloc (sizeof (interval_t) * index->read_len * 2);
        ts->I_B =
            (interval_t *) malloc (sizeof (interval_t) * index->read_len * 2);
        assert (ts->I_B != NULL);
        assert (ts->I_F != NULL);
        assert (ts->I_fw != NULL);
        assert (ts->I_bw != NULL);
        ts->I_F[0].uu = index->read_len * 2;
        ts->I_B[0].uu = index->read_len * 2;
        ts->I_B[1].uu = 0;
        ts->I_F[1].uu = 0;

        INDEX_TYPE_T i;

#ifndef DEBUG1
#ifdef __OPENMP__
        #pragma omp for schedule(dynamic, 1) reduction(+: overlaps) reduction(+: overlap_len)
#endif
#endif
        // parallel for each read
        for (i = 0; i < index->num_reads; i++)
        {
            int j;
            int onum;
            int olen;
            node_t *node;
            edge_t *new_edge;
            int overlaps_B;
            int overlaps_F;
            double ratio;
            path_q_t *path_q;
            path_t *path;
            ts->pos_fw = 0;
            ts->pos_bw = 0;
            ts->size_I_F = 0;
            ts->size_I_B = 0;
            onum = 0;
            olen = 0;
                    
            // find overlapping candidates
            find_intervals (i, index, ts, &overlaps_B, &overlaps_F);
           
            // remove transitive edges
            ts->level = 0;
            if (overlaps_B != 0)
                extract_bw_edges (graph, i, index, ts);

            ts->level = 0;
            if (overlaps_F != 0)
                extract_fw_edges (graph, i, index, ts);
            
            // add node and edges           
            node = &(graph->nodes[i]);
            node->max_B = node->num_B = ts->size_I_B;
            node->max_F = node->num_F = ts->size_I_F;
            node->status = DEFAULT;
            node->path_q = NULL;
            node->edges = (edge_t *) malloc ((ts->size_I_B + ts->size_I_F) * sizeof (edge_t));
            assert (node->edges);
            olen = 0;
            onum = 0;

            // marked for bridge removal
            ratio = max(overlaps_B, overlaps_F) == 0 ?
                    0 : ((double)abs(overlaps_B - overlaps_F))/
                          max(overlaps_B, overlaps_F);

            if (ratio > 0.6)
            {
                node->marked = SUSPECT;               
            }
            else
            {
                node->marked = 0;
            }

            // classify the nodes
            if (node->num_B != 0 || node->num_F != 0)
            { 
                if (node->num_B == 1)
                {
                    node->B_type = ONE;
                }
                else if (node->num_B > 1)
                {
                    node->B_type = BRANCH;
                }
                else
                {
                    node->B_type = ENDING;
                }

                if (node->num_F == 1)
                {
                    node->F_type = ONE;
                }
                else if (node->num_F > 1)
                {
                    node->F_type = BRANCH;
                }
                else
                {
                    node->F_type = ENDING;
                }
            }
            else
            {
                node->B_type = node->F_type = ENDING;
                node->status = REMOVED;
            }

#ifdef DEBUG1
            fprintf (fp, "read %d: <%d %d> ", i, overlaps_B, overlaps_F);
            fprintf (fp1, "%d %d %d %lf %lf\n", i, overlaps_B, overlaps_F);
#endif

            // malloc edges and paths
            if ((node->B_type == BRANCH && node->F_type != ENDING) ||
                (node->F_type == BRANCH && node->B_type != ENDING))
            {
                path_q = node->path_q = (path_q_t *) malloc (sizeof (path_q_t));
                assert (path_q != NULL);
                path_q->max_B_endlen = 0;
                path_q->max_F_endlen = 0;
                path_q->max_B_len = 0;
                path_q->max_F_len = 0;
                path_q->path = (path_t *) malloc (sizeof (path_t) *
                                                  (node->max_B + node->max_F));
                path_q->max_F_trans = ts->I_F[1].uu;
                path_q->max_B_trans = ts->I_B[1].uu;
                assert (path_q->path != NULL);

                for (j = 0; j < node->max_B; j++)
                {
                    new_edge = &(node->edges[j]);
                    new_edge->dest = ts->I_B[j].ll;
                    new_edge->overlap_len = ts->I_B[j].len;
                    new_edge->type = ts->I_B[j].type;
                    onum++;
                    olen += new_edge->overlap_len;
#ifdef DEBUG1
                    fprintf (fp, "(%d, ", new_edge->dest);
                    fprintf (fp, "%d, ", new_edge->overlap_len);
                    if (new_edge->type == BF)
                    {
                        fprintf (fp, "BF); ");
                    }
                    else
                    {
                        fprintf (fp, "BB); ");
                    }
#endif
                    path = &(path_q->path[j]);
                    path->end = new_edge->dest;
                    path->hop = new_edge->dest;
                    path->len = 1;
                    path->end_len = graph->read_len - new_edge->overlap_len;
                    path->bubble_len = 0;
                    path->type = NONE;
                    path->direction = (new_edge->type == BB ? 1 : 0);
                    path->trans = ts->I_B[j].remain_size;
                }
                for (j = 0; j < node->max_F; j++)
                {
                    new_edge = &(node->edges[j + node->max_B]);
                    new_edge->dest = ts->I_F[j].ll;
                    new_edge->overlap_len = ts->I_F[j].len;
                    new_edge->type = ts->I_F[j].type;
                    onum++;
                    olen += new_edge->overlap_len;
#ifdef DEBUG1
                    fprintf (fp, "(%d, ", new_edge->dest);
                    fprintf (fp, "%d, ", new_edge->overlap_len);
                    if (new_edge->type == FF)
                    {
                        fprintf (fp, "FF); ");
                    }
                    else
                    {
                        fprintf (fp, "FB); ");
                    }
#endif
                    path = &(path_q->path[j + node->max_B]);
                    path->end = new_edge->dest;
                    path->hop = new_edge->dest;
                    path->len = 1;
                    path->end_len = graph->read_len - new_edge->overlap_len;
                    path->bubble_len = 0;
                    path->type = NONE;
                    path->direction = (new_edge->type == FB ? 1 : 0);
                    path->trans = ts->I_F[j].remain_size;
                }
            }
            else
            {
                for (j = 0; j < node->max_B; j++)
                {
                    new_edge = &(node->edges[j]);
                    new_edge->dest = ts->I_B[j].ll;
                    new_edge->overlap_len = ts->I_B[j].len;
                    new_edge->type = ts->I_B[j].type;
                    onum++;
                    olen += new_edge->overlap_len;
#ifdef DEBUG1
                    fprintf (fp, "(%d, ", new_edge->dest);
                    fprintf (fp, "%d, ", new_edge->overlap_len);
                    if (new_edge->type == BF)
                    {
                        fprintf (fp, "BF); ");
                    }
                    else
                    {
                        fprintf (fp, "BB); ");
                    }
#endif
                }
                for (j = 0; j < node->max_F; j++)
                {
                    new_edge = &(node->edges[j + node->max_B]);
                    new_edge->dest = ts->I_F[j].ll;
                    new_edge->overlap_len = ts->I_F[j].len;
                    new_edge->type = ts->I_F[j].type;
                    onum++;
                    olen += new_edge->overlap_len;
#ifdef DEBUG1
                    fprintf (fp, "(%d, ", new_edge->dest);
                    fprintf (fp, "%d, ", new_edge->overlap_len);
                    if (new_edge->type == FF)
                    {
                        fprintf (fp, "FF); ");
                    }
                    else
                    {
                        fprintf (fp, "FB); ");
                    }
#endif
                }
            }
        
           
            overlaps += onum;
            overlap_len += (double) olen;
#ifdef DEBUG1
            fprintf (fp, "end\n");
#endif
        }
        free (ts->I_bw);
        free (ts->I_fw);
        free (ts->I_B);
        free (ts->I_F);
        free (ts);

#ifndef DEBUG1
#ifdef __OPENMP__
    }
#endif
#endif

#ifdef __TIMING__
    gettimeofday (&tv2, NULL);
    time_pass =
        (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec -
                                              tv1.tv_usec) / 1000.0;
    DPRINTF (3, "  takes %.3lf ms\n", time_pass);
#endif

#ifdef DEBUG1
    fclose (fp);
    fclose (fp1);
#endif

    graph->average_overlap = overlap_len / overlaps;
    DPRINTF (3, "  average overlap length: %lf\n", graph->average_overlap);
    return graph;    
}


void
destroy_graph (string_graph_t * graph)
{
}
