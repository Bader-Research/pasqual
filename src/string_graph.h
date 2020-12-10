#ifndef __OVERLAP_GRAPH_H__
#define __OVERLAP_GRAPH_H__


#include "common.h"
#include "index.h"


/* -------------------------------------------------------------------------- */
/* enums                                                                      */

// edge type
#define BF      (char)0
#define FB      (char)1
#define FF      (char)2
#define BB      (char)3
#define INVALID (char)4

// node type
#define ONE     (char)1
#define BRANCH  (char)2
#define HALF    (char)3
#define ENDING  (char)4
#define SELF    (char)5
#define TENDING (char)6

// node status
#define DEFAULT   (char)0
#define VISITED   ((char)1)<<1
#define CHANGED   ((char)1)<<2
#define REMOVED   ((char)1)<<3
#define CHANGED2  ((char)1)<<4

// path type
#define TIP      (char)1
#define BUBBLE   (char)2
#define NONE     (char)0

// marked type
#define MARKED   ((char)1)<<1
#define SUSPECT  ((char)1)<<2

/* -------------------------------------------------------------------------- */
/* type of structures                                                         */

// suffix array interval
typedef struct _interval_t
{
    INDEX_TYPE_T ll;
    INDEX_TYPE_T uu;
    INDEX_TYPE_T remain_size;
    u16 len;
    char type;
} interval_t;

// path
typedef struct _path_t
{
    INDEX_TYPE_T end;
    INDEX_TYPE_T hop;
    u16 len;
    u16 end_len;
    u16 bubble_len;
    int trans;
    char updated;
    char direction;             // 1: F, 0: B
    char type;                  // 1: TIP,  0: BUBBLE
} path_t;

// path list
typedef struct _path_q_t
{
    path_t *path;
    u16 max_F_endlen;
    u16 max_B_endlen;
    u16 max_B_len;
    u16 max_F_len;
    int max_F_trans;
    int max_B_trans;
} path_q_t;

// edge
typedef struct _edge_t
{
    INDEX_TYPE_T dest;
    u16 overlap_len;
    char type;
} edge_t;

// node
typedef struct _node_t
{
    edge_t *edges;
    path_q_t *path_q;
    u16 num_F;
    u16 num_B;
    u16 max_F;
    u16 max_B;
    char B_type;
    char F_type;
    char status;
    char marked;
} node_t;

// string graph
typedef struct _string_graph_t
{
    u16 read_len;
    u16 min_overlap;
    double average_overlap;
    INDEX_TYPE_T num_reads;
    node_t *nodes;
} string_graph_t;

// resource of thread
typedef struct _thread_resource_t
{
    interval_t *I_B;
    interval_t *I_F;
    interval_t *I_fw;
    interval_t *I_bw;
    int pos_fw;
    int pos_bw;
    int size_fw;
    int size_bw;
    int level;
    u16 size_I_F;
    u16 size_I_B;
    u16 max_fw;
    u16 max_bw;
    int max_trans;
} thread_resource_t;


/* -------------------------------------------------------------------------- */
/* definition of functions                                                    */

// build string graph
string_graph_t *construct_string_graph (index_t * index);

// destroy string graph
void destroy_graph (string_graph_t * graph);


#endif /* __OVERLAP_GRAPH_H__ */
