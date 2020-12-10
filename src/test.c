#include <stdio.h>
#include <stdlib.h>
#ifdef __TIMING__
#include <time.h>
#include <sys/time.h>
#endif
#include <getopt.h>

#include "index.h"
#include "string_graph.h"
#include "common.h"
#include "assemble.h"
#include "bloom.h"

/* -------------------------------------------------------------------------- */
/* Open issues or limits:                                                     */
/* 1. Length of reads is shorter than 16 bits (65535)                         */
/* 2. Number of overlaps of each read is less than 16 bits (65535)            */
/* -------------------------------------------------------------------------- */

const struct option long_options[] = {
    {"help", 0, NULL, 'h'},
    {"version", 0, NULL, 'v'},
    {"singlereads", 1, NULL, 11},
#if 0
    {"paired", 1, NULL, 12},
    {"mate", 1, NULL, 14},
#endif
    {"overlap", 1, NULL, 't'},
    {"readlen", 1, NULL, 'l'},
    {"maxtiplen", 1, NULL, 16},
    {"mincontiglen", 1, NULL, 17},
    {"bloomlen", 1, NULL, 18},
    {"nhashfunc", 1, NULL, 19},
    {"numiters", 1, NULL, 20},    
    {"bubblelen", 1, NULL, 21},
    {NULL, 0, NULL, 0}
};

const char *const short_options = "-:hvt:l:";
const char *version_info = "1.0.0";


// print usage
static void
usage (char *call)
{
    DPRINTF (1, "Usage: %s [OPTIONS]\n", call);
    DPRINTF (1, "Options:\n");
    DPRINTF (1, "\t-h or --help        Display this information\n");
    DPRINTF (1, "\t-v or --version     Display version information\n");
    DPRINTF (1, "\t-l or --readlen     Reads length\n");
    DPRINTF (1, "\t-t or --overlap     Minimum overlaps length\n");
    DPRINTF (1,
             "\t--mincontiglen      Minimum contigs length, default = 100\n");
    DPRINTF (1,
             "\t--maxtiplen         Maximum tip length, default = 15\n");
    DPRINTF (1,
             "\t--bubblelen         Bubbles length, default = 10\n");
    DPRINTF (1,
             "\t--numiters          Number of iterations for tips removal, default = maxtiplen\n");
    DPRINTF (1,
             "\t--bloomlen          Length of bloom filter array, default = 36\n");
    DPRINTF (1,
             "\t--nhashfunc         Number of hash functions used for bloom filter, default = 2\n");
    DPRINTF (1, "\t--singlereads       Single reads input in FASTA format\n");
#if 0
    DPRINTF (1,
             "\t--paired            Paired-end reads input (2 files) in FASTA format\n");
    DPRINTF (1,
             "\t--mate              Mate-pair reads input (2 files) in FASTA format\n");
#endif
}


static void
version (char *call)
{
    DPRINTF (1, "%s version %s\n", call, version_info);
}


// main function
int
main (int argc, char **argv)
{
#ifdef __TIMING__
    struct timeval tv1, tv2;
    double time_pass;
#endif

    index_t *index;
    string_graph_t *graph;

    // read length
    int read_len;

    // minimum overlap length
    int min_overlap;

    // minimum length of contigs
    int min_contiglen;

    // maximum tip length
    int max_tiplen;

    // bubble length
    int bubblelen;
    
    // number of contigs
    INDEX_TYPE_T num_contigs;

    // total length of contigs
    INDEX_TYPE_T len_contigs;

    // length of bloom filter array
    int bloom_m;

    // number of hash functions of bloom filter
    int n_hashfunc;

    // number iterations for tips removal
    int num_iters;

    int i;
    int c = 0;
    int pflag;
    int mflag;
    char *single_file;
    char *pair_file1;
    char *pair_file2;
    char *mate_file1;
    char *mate_file2;
    char *call;

    // default parameters
    // read len is less than 16bits, minimum overlap is longer than SEQ_DEPTH    
    read_len = -1;
    min_overlap = -1;
    num_iters = 0;
    max_tiplen = 0;
    min_contiglen = 0;
    bloom_m = 36;
    n_hashfunc = 2;
    bubblelen = 10;
    single_file = NULL;
    pair_file1 = pair_file2 = NULL;
    mate_file1 = mate_file2 = NULL;

    call = argv[0];
    while (strstr (call, "/") != NULL)
    {
        call = strstr (call, "/");
        call++;
    }    
    pflag = 0;
    mflag = 0;
    /* parse arguments */
    while ((c =
            getopt_long (argc, argv, short_options, long_options,
                         NULL)) != -1)
    {
        switch (c)
        {
        case 'h':
            usage (call);
            return 0;
        case 'v':
            version (call);
            return 0;
        case 'l':
            read_len = atoi (optarg);
            if (read_len <= 0)
            {
                ERROR_EXIT ("%s: reads length must be larger than 0.\n",
                            call);
            }
            break;
        case 't':
            min_overlap = atoi (optarg);
            if (min_overlap <= 0)
            {
                ERROR_EXIT
                    ("%s: minimum overlaps length must be larger than 0.\n",
                     call);
            }
            break;
        case 1:
            if (pflag == 1)
            {
                pair_file2 =
                    (char *) malloc ((strlen (optarg) + 1) * sizeof (char));
                assert (pair_file2 != NULL);
                strcpy (pair_file2, optarg);
                pflag = 0;
            }
            else if (mflag == 1)
            {
                mate_file2 =
                    (char *) malloc ((strlen (optarg) + 1) * sizeof (char));
                assert (mate_file2 != NULL);
                strcpy (mate_file2, optarg);
                mflag = 0;
            }
            else
            {
                DPRINTF (1, "%s: non-option ARGV-elements: ", call);
                while (optind <= argc)
                {
                    DPRINTF (1, "%s ", argv[optind - 1]);
                    optind++;
                }
                DPRINTF (1, "\n");
                DPRINTF (1, "Try %s --help for more information.\n", call);
                return -1;
            }
            break;
        case 11:
            single_file =
                (char *) malloc ((strlen (optarg) + 1) * sizeof (char));
            assert (single_file != NULL);
            strcpy (single_file, optarg);
            break;
        case 12:
            pair_file1 =
                (char *) malloc ((strlen (optarg) + 1) * sizeof (char));
            assert (pair_file1 != NULL);
            strcpy (pair_file1, optarg);
            pflag = 1;
            break;
        case 14:
            mate_file1 =
                (char *) malloc ((strlen (optarg) + 1) * sizeof (char));
            assert (mate_file1 != NULL);
            strcpy (mate_file1, optarg);
            mflag = 1;
            break;
        case 16:
            max_tiplen = atoi (optarg);
            break;
        case 17:
            min_contiglen = atoi (optarg);
            break;
        case 18:
            bloom_m = atoi (optarg);
            break;
        case 19:
            n_hashfunc = atoi (optarg);
            break;
        case 20:
            num_iters = atoi (optarg);
            break;
        case 21:
            bubblelen = atoi (optarg);
            break;            
        case ':':
            DPRINTF (1, "%s: option '%s' requires an argument.\n", call,
                     argv[optind - 1]);
            DPRINTF (1, "Try %s --help for more information.\n", call);
            return -1;
        case '?':
            if (argv[optind - 1][0] != '-')
            {
                DPRINTF (1, "%s: unknown option `%s'.\n", call,
                         argv[optind]);
            }
            else
            {
                DPRINTF (1, "%s: unknown option `%s'.\n", call,
                         argv[optind - 1]);
            }
            DPRINTF (1, "Try %s --help for more information.\n", call);
            return -1;
        default:
            usage (call);
            return -1;
        }
    }

    // check parameters
    if (read_len == -1)
    {
        ERROR_EXIT ("missing required argument --readlen(-l).\n");
    }
    if (min_overlap == -1)
    {
        ERROR_EXIT ("missing required argument --overlap(-t).\n");
    }
    if (pair_file1 != NULL && pair_file2 == NULL)
    {
        ERROR_EXIT ("paired-end reads must be two files.\n");
    }
    if (mate_file1 != NULL && mate_file2 == NULL)
    {
        ERROR_EXIT ("mate-pair reads must be two files.\n");
    }
    if (single_file == NULL && pair_file1 == NULL && mate_file1 == NULL)
    {
        ERROR_EXIT ("missing reads input.\n");
    }
    if (read_len <= 0 || read_len > 65535)
    {
        ERROR_EXIT ("invalid read length (%d), valid range (0, 65535].\n",
                    read_len);
    }
    if (min_overlap <= SEQ_DEPTH || min_overlap >= read_len)
    {
        ERROR_EXIT
            ("invalid minimum overlaps length(%d), valid range (%d, %d).\n",
             min_overlap, SEQ_DEPTH, read_len);
    }
    if (n_hashfunc < 1 || n_hashfunc > MAX_FUNC)
    {
        ERROR_EXIT
            ("invalid number of hash functions (%d), valid range [1, MAX_FUNC].\n",
             n_hashfunc);
    }
    if (bloom_m < 1)
    {
        ERROR_EXIT ("length of bloom filter array must be larger than 1.\n");
    }
    if (max_tiplen < 0)
    {
        ERROR_EXIT ("maximum tip length must be larger than 0.\n");
    }
    if (bubblelen <= 0)
    {
        ERROR_EXIT ("bubble length must be larger than 0.\n");
    }
    if (num_iters < 0 || num_iters > max_tiplen)
    {
        ERROR_EXIT
            ("number interations must be larger than 0 and smaller than maxtiplen.\n");
    }

    if (num_iters == 0)
    {
        num_iters = 2*read_len/5;
    }
    if (max_tiplen == 0)
    {
        max_tiplen = 2*read_len;
    }
    if (min_contiglen == 0)
    {
        min_contiglen = 2*read_len;
    }
    
    version ("PASQUAL");
    
    DPRINTF (1, "  input reads:\t%s\n", single_file);
    DPRINTF (1, "  output:\t%s\n", "contig.fa"); 
    DPRINTF (1, "  read length:\t%d\n", read_len);
    DPRINTF (1, "  min overlap:\t%d\n", min_overlap);    
    DPRINTF (1, "  min contig:\t%d\n", min_contiglen);
    DPRINTF (1, "  max tip:\t%d\n", max_tiplen);
    DPRINTF (1, "  bubble:\t%d\n", bubblelen);
        
#ifdef __TIMING__
    gettimeofday (&tv1, NULL);
#endif

    // STEP 1: build compressed index
    index = construct_index (single_file,
                             pair_file1, pair_file2,
                             mate_file1, mate_file2,
                             min_overlap, read_len, bloom_m, n_hashfunc);

    // STEP 2: build string graph
    graph = construct_string_graph (index);

    free (index->OCC_sample);
    free (index->OCC_sample_bw);
    for (i = 0; i < 5; i++)
    {
        free (index->OCC_bitmap[i]);
        free (index->OCC_bitmap_bw[i]);
    }

    free (index->reads_index);
    free (index->reads_index_bw);
    free (index->interval_depth);
    free (index->interval_depth_bw);
    free (index->interval_end);
    free (index->interval_end_bw);
    free (index->interval_start);
    free (index->interval_start_bw);

    // STEP 3: list contigs
    assemble (graph, index->reads, &num_contigs, &len_contigs, min_contiglen,
              max_tiplen, bubblelen, num_iters);

    DPRINTF (1, "Summary:\n");
#if defined(__INDEX_U32__)
    DPRINTF (1, "  number of contigs:\t%d\n", num_contigs);
    DPRINTF (1, "  total contig length:\t%d\n", len_contigs);
#elif defined(__INDEX_U64__)
    DPRINTF (1, "  number of contigs:\t%lld\n", num_contigs);
    DPRINTF (1, "  total contig length:\t%lld\n", len_contigs);
#else
    DPRINTF (1, "  number of contigs:\t%d\n", num_contigs);
    DPRINTF (1, "  total contig length:\t%d\n", len_contigs);
#endif

    free (index->reads);
    free (index);
    destroy_graph (graph);
    free (graph);

#ifdef __TIMING__
    gettimeofday (&tv2, NULL);
    time_pass =
        (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec -
                                              tv1.tv_usec) / 1000.0;
    DPRINTF (1,
             "  total time:\t%.3lf ms\n", time_pass);
    DPRINTF (1,
             "  exclude I/O:\t%.3lf ms\n",
             time_pass - time_importreads);
#endif

    if (single_file != NULL)
        free (single_file);
    if (pair_file1 != NULL)
        free (pair_file1);
    if (pair_file2 != NULL)
        free (pair_file2);
    if (mate_file1 != NULL)
        free (mate_file1);
    if (mate_file2 != NULL)
        free (mate_file2);

    return 0;
}
