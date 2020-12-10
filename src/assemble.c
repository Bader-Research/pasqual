#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "string_graph.h"
#ifdef __OPENMP__
#include <omp.h>
#endif
#ifdef __TIMING__
#include <time.h>
#include <sys/time.h>
#endif


// TODO: clean the memory

inline void
destroy_node (node_t * node)
{
    if (node->path_q->path != NULL)
    {
        free (node->path_q->path);
        node->path_q->path = NULL;
    }
    if (node->edges != NULL)
    {
        free (node->edges);
        free (node->path_q);
        node->edges = NULL;
        node->path_q = NULL;
    }
}


#ifdef DEBUG2
static void
check_missing_edges (string_graph_t * graph)
{
    INDEX_TYPE_T i;
    int num_miss = 0;
    INDEX_TYPE_T nnnn;
    nnnn = 0;

    for (i = 0; i < graph->num_reads; i++)
    {
        node_t *node = &(graph->nodes[i]);
        edge_t *d_edge;
        node_t *d_node;
        edge_t *edge;
        int flag;
        int k;
        int j;

        if ((node->status & REMOVED) != 0)
        {
            continue;
        }

        // B-edge side
        for (j = 0; j < node->max_B; j++)
        {
            edge = &(node->edges[j]);
            if (edge->type == INVALID)
                continue;
            nnnn++;
            d_node = &(graph->nodes[edge->dest]);
            flag = 0;
            // has a BF edge
            if (edge->type == BF)
            {
                // check whether d_node has a FB edge
                for (k = 0; k < d_node->max_F; k++)
                {
                    d_edge = &(d_node->edges[k + d_node->max_B]);
                    if (d_edge->type == FB &&
                        d_edge->dest == i &&
                        d_edge->overlap_len == edge->overlap_len)
                    {
                        flag = 1;
                        break;
                    }
                }
                // miss
                if (flag == 0)
                {
                    num_miss++;
#if defined(__INDEX_U32__)
                    DPRINTF (1, "%d %d FB %d\n", edge->dest, i,
                             edge->overlap_len);
#elif defined(__INDEX_U64__)
                    DPRINTF (1, "%lld %lld FB %d\n", edge->dest, i,
                             edge->overlap_len);
#else
                    DPRINTF (1, "%d %d FB %d\n", edge->dest, i,
                             edge->overlap_len);
#endif
                }
            }
            // has a BB edge
            else
            {
                // check whether d_node has a BB edge
                for (k = 0; k < d_node->max_B; k++)
                {
                    d_edge = &(d_node->edges[k]);
                    if (d_edge->type == BB &&
                        d_edge->dest == i &&
                        d_edge->overlap_len == edge->overlap_len)
                    {
                        flag = 1;
                        break;
                    }
                }
                // miss
                if (flag == 0)
                {
                    num_miss++;
#if defined(__INDEX_U32__)
                    DPRINTF (1, "%d %d BB %d\n", edge->dest, i,
                             edge->overlap_len);
#elif defined(__INDEX_U64__)
                    DPRINTF (1, "%lld %lld BB %d\n", edge->dest, i,
                             edge->overlap_len);
#else
                    DPRINTF (1, "%d %d BB %d\n", edge->dest, i,
                             edge->overlap_len);
#endif
                }
            }

        }

        // F-edge side
        for (j = 0; j < node->max_F; j++)
        {
            edge = &(node->edges[j + node->max_B]);
            if (edge->type == INVALID)
                continue;
            nnnn++;
            d_node = &(graph->nodes[edge->dest]);
            flag = 0;
            // has a FF edge
            if (edge->type == FF)
            {
                // check whether d_node has a FF edge
                for (k = 0; k < d_node->max_F; k++)
                {
                    d_edge = &(d_node->edges[k + d_node->max_B]);
                    if (d_edge->type == FF &&
                        d_edge->dest == i &&
                        d_edge->overlap_len == edge->overlap_len)
                    {
                        flag = 1;
                        break;
                    }
                }
                // miss
                if (flag == 0)
                {
                    num_miss++;
#if defined(__INDEX_U32__)
                    DPRINTF (1, "%d %d FF %d\n", edge->dest, i,
                             edge->overlap_len);
#elif defined(__INDEX_U64__)
                    DPRINTF (1, "%lld %lld FF %d\n", edge->dest, i,
                             edge->overlap_len);
#else
                    DPRINTF (1, "%d %d FF %d\n", edge->dest, i,
                             edge->overlap_len);
#endif
                }
            }
            // has a FB edge
            else
            {
                // check whether d_node has a BF edge
                for (k = 0; k < d_node->max_B; k++)
                {
                    d_edge = &(d_node->edges[k]);
                    if (d_edge->type == BF &&
                        d_edge->dest == i &&
                        d_edge->overlap_len == edge->overlap_len)
                    {
                        flag = 1;
                        break;
                    }
                }
                // miss
                if (flag == 0)
                {
                    num_miss++;
#if defined(__INDEX_U32__)
                    DPRINTF (1, "%d %d BF %d\n", edge->dest, i,
                             edge->overlap_len);
#elif defined(__INDEX_U64__)
                    DPRINTF (1, "%lld %lld BF %d\n", edge->dest, i,
                             edge->overlap_len);
#else
                    DPRINTF (1, "%d %d BF %d\n", edge->dest, i,
                             edge->overlap_len);
#endif
                }
            }
        }
    }
    DPRINTF (1, "miss %d edges %d ", num_miss, nnnn);
}
#endif


#ifdef DEBUG2
// remove bubble path
static void
remove_path0 (INDEX_TYPE_T head, int direction, string_graph_t * graph)
{
    node_t *node;
    int i;
    edge_t *edge;
    INDEX_TYPE_T end;

    end = head;
    node = &(graph->nodes[end]);

    // traverse until an ENDING node
    while (node->B_type == ONE && node->F_type == ONE)
    {
        node->status = REMOVED;
        if (direction == 1)
        {
            for (i = 0; i < node->max_F; i++)
            {
                edge = &(node->edges[i + node->max_B]);
                if (edge->type != INVALID)
                {
                    end = edge->dest;
                    direction = (edge->type == FB ? 1 : 0);
                    break;
                }
            }
        }
        else
        {
            for (i = 0; i < node->max_B; i++)
            {
                edge = &(node->edges[i]);
                if (edge->type != INVALID)
                {
                    end = edge->dest;
                    direction = (edge->type == BB ? 1 : 0);
                    break;
                }
            }
        }
        node = &(graph->nodes[end]);
    }
}


// remove junction path
static void
remove_path1 (INDEX_TYPE_T head, int direction, string_graph_t * graph)
{
    node_t *node;
    int i;
    edge_t *edge;
    INDEX_TYPE_T end;

    end = head;
    node = &(graph->nodes[end]);

    // traverse until an ENDING node
    while (node->B_type == ONE && node->F_type == ONE)
    {
        node->status = REMOVED;
        if (direction == 1)
        {
            for (i = 0; i < node->max_F; i++)
            {
                edge = &(node->edges[i + node->max_B]);
                if (edge->type != INVALID)
                {
                    end = edge->dest;
                    direction = (edge->type == FB ? 1 : 0);
                    break;
                }
            }
        }
        else
        {
            for (i = 0; i < node->max_B; i++)
            {
                edge = &(node->edges[i]);
                if (edge->type != INVALID)
                {
                    end = edge->dest;
                    direction = (edge->type == BB ? 1 : 0);
                    break;
                }
            }
        }
        node = &(graph->nodes[end]);
    }
    node->marked |= MARKED;
}


// remove tip path
static void
remove_path3 (INDEX_TYPE_T head, int direction, string_graph_t * graph)
{
    node_t *node;
    node_t *d_node;
    int i;
    edge_t *edge;
    INDEX_TYPE_T end;

    end = head;
    node = &(graph->nodes[end]);
    node->status = REMOVED;

    // traverse until an ENDING node
    while (node->B_type == ONE && node->F_type == ONE)
    {
        if (direction == 1)
        {
            for (i = 0; i < node->max_F; i++)
            {
                edge = &(node->edges[i + node->max_B]);
                if (edge->type != INVALID)
                {
                    end = edge->dest;
                    direction = (edge->type == FB ? 1 : 0);
                    break;
                }
            }
        }
        else
        {
            for (i = 0; i < node->max_B; i++)
            {
                edge = &(node->edges[i]);
                if (edge->type != INVALID)
                {
                    end = edge->dest;

                    direction = (edge->type == BB ? 1 : 0);
                    break;
                }
            }
        }
        node = &(graph->nodes[end]);
        node->status = REMOVED;
    }

    // handle dead end
    if (direction == 0 && node->F_type == BRANCH)
    {
        for (i = 0; i < node->max_F; i++)
        {
            edge = &(node->edges[i + node->max_B]);
            if (edge->type != INVALID)
            {
                d_node = &(graph->nodes[edge->dest]);
                if (d_node->B_type == ENDING ||
                     d_node->F_type == ENDING ||
                     (d_node->B_type == ONE && d_node->F_type == ONE))
                {
                    d_node->marked |= MARKED;
                }
            }
        }
    }
    else if (direction == 1 && node->B_type == BRANCH)
    {
        for (i = 0; i < node->max_B; i++)
        {
            edge = &(node->edges[i]);
            if (edge->type != INVALID)
            {
                d_node = &(graph->nodes[edge->dest]);
                if (d_node->B_type == ENDING ||
                     d_node->F_type == ENDING ||
                     (d_node->B_type == ONE && d_node->F_type == ONE))
                {
                    d_node->marked |= MARKED;
                }
            }
        }
    }
}

#else

// remove junction path
static void
remove_path2 (INDEX_TYPE_T head, path_t * path, string_graph_t * graph)
{
    node_t *node;
    if (head != path->end)
    {
        node = &(graph->nodes[head]);
        node->status = REMOVED;
        node = &(graph->nodes[path->hop]);
        node->status = REMOVED;
    }        
    node = &(graph->nodes[path->end]);
    node->marked |= MARKED;
}


// remove tip path
static void
remove_path4 (INDEX_TYPE_T head, path_t * path,
              int direction, string_graph_t * graph)
{
    node_t *node;
    node_t *head_node;
    node_t *hop_node;
    node_t *d_node;
    int i;
    edge_t *edge;

    node = &(graph->nodes[path->end]);
    
    // handle dead end
    if (direction == 0 && node->F_type == BRANCH)
    {
        for (i = 0; i < node->max_F; i++)
        {
            edge = &(node->edges[i + node->max_B]);
            if (edge->type != INVALID)
            {
                d_node = &(graph->nodes[edge->dest]);
                if ((node->status & REMOVED) == 0 &&
                    (d_node->B_type == ENDING ||
                     d_node->F_type == ENDING ||
                     (d_node->B_type == ONE && d_node->F_type == ONE)))
                {
                    d_node->marked |= MARKED;
                }
            }
        }
    }
    else if (direction == 1 && node->B_type == BRANCH)
    {
        for (i = 0; i < node->max_B; i++)
        {
            edge = &(node->edges[i]);
            if (edge->type != INVALID)
            {
                d_node = &(graph->nodes[edge->dest]);
                if ((node->status & REMOVED) == 0 &&
                    (d_node->B_type == ENDING ||
                     d_node->F_type == ENDING ||
                     (d_node->B_type == ONE && d_node->F_type == ONE)))
                {
                    d_node->marked |= MARKED;
                }
            }
        }
    }
    
    head_node = &(graph->nodes[head]);
    head_node->status = REMOVED;
    hop_node = &(graph->nodes[path->hop]);
    hop_node->status = REMOVED;    
    node->status = REMOVED;
}
#endif


// extend path
static void
extend_path (path_t * path, string_graph_t * graph)
{
    node_t *node;
    int direction;
    int i;
    edge_t *n_edge;
    INDEX_TYPE_T end;
    INDEX_TYPE_T hop;
    INDEX_TYPE_T len;
    INDEX_TYPE_T end_len;
    INDEX_TYPE_T bubble_len;
    int read_len;

    if (path->type == TIP)
        return;

    end = path->end;
    hop = path->hop;
    len = path->len;
    end_len = path->end_len;
    direction = path->direction;
    node = &(graph->nodes[end]);
    read_len = graph->read_len;
    n_edge = NULL;

    while (node->B_type == ONE && node->F_type == ONE)
    {
        hop = end;
        if (direction == 1)
        {
            for (i = 0; i < node->max_F; i++)
            {
                n_edge = &(node->edges[i + node->max_B]);
                if (n_edge->type != INVALID)
                {
                    end = n_edge->dest;
                    end_len += (read_len - n_edge->overlap_len);
                    direction = (n_edge->type == FB ? 1 : 0);
                    break;
                }
            }
        }
        else
        {
            for (i = 0; i < node->max_B; i++)
            {
                n_edge = &(node->edges[i]);
                if (n_edge->type != INVALID)
                {
                    end = n_edge->dest;
                    end_len += (read_len - n_edge->overlap_len);
                    direction = (n_edge->type == BB ? 1 : 0);
                    break;
                }
            }
        }
        node = &(graph->nodes[end]);
        len++;
    }
    path->end = end;
    path->hop = hop;
    len = len > 65535 ? 63535 : len;
    path->len = len;
    path->updated = 1;
    path->direction = direction;
    if (n_edge != NULL)
    {
        bubble_len = end_len - (read_len - n_edge->overlap_len);
        bubble_len = bubble_len > 65535 ? 63535 : bubble_len;
        path->bubble_len = bubble_len;
    }
    end_len = end_len > 65535 ? 63535 : end_len;
    path->end_len = end_len;

    /* identify path type */
    if (node->B_type == ENDING || node->F_type == ENDING)
    {
        path->type = TIP;
    }
    else if ((node->B_type == BRANCH && direction == 1) ||
             (node->F_type == BRANCH && direction == 0))
    {
        path->type = BUBBLE;
    }
    else
    {
        path->type = NONE;
    }
}


// remove bridge
static void
remove_bridge (string_graph_t * graph)
{
#ifdef DEBUG2
    FILE *fp;
    fp = fopen ("filter.dat", "w+");
#endif

    INDEX_TYPE_T i;
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {
        int k;
        node_t *node;
        edge_t *edge;
        path_q_t *path_q;
        path_t *path;
        int flag;
#ifdef DEBUG2
        int direction;
#endif
        node = &(graph->nodes[i]);
        path_q = node->path_q;
        flag = 0;
        
        // find bridge edges and nodes
        if ( (node->status & REMOVED) == 0 &&
             (node->marked & SUSPECT) != 0 &&
             ((node->B_type == BRANCH && node->F_type != ENDING) ||
              (node->F_type == BRANCH && node->B_type != ENDING)) )
        {
            for (k = 0; k < node->max_F; k++)
            {
                edge = &(node->edges[k + node->max_B]);
                path = &(path_q->path[k + node->max_B]);
                if (edge->type != INVALID &&
                    (path->type == TIP ||
                     path->len >= 3))
                {                  
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    if (edge->type != INVALID &&
                        (path->type == TIP ||
                         path->len >= 3))
                    {                  
                        flag = 1;
                        break;
                    }
                }
            }

            // remove bridge
            if (flag == 0)
            {
#ifdef DEBUG2            
                fprintf (fp, "%d\n", i);
#endif
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    if (edge->type != INVALID)
                    {
#ifdef DEBUG2
                        direction = (edge->type == BB ? 1 : 0);
                        remove_path1 (edge->dest, direction, graph);
#else
                        remove_path2 (edge->dest, path, graph);
#endif
                    }
                }
                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    path = &(path_q->path[k + node->max_B]);
                    if (edge->type != INVALID)
                    {
#ifdef DEBUG2
                        direction = (edge->type == FB ? 1 : 0);
                        remove_path1 (edge->dest, direction, graph);
#else
                        remove_path2 (edge->dest, path, graph);
#endif
                    }
                }
                node->status = REMOVED;
            }
        }
    }

    // remove from the other end
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {
        int k;
        node_t *node;
        edge_t *edge;
        node_t *d_node;
        
        node = &(graph->nodes[i]);
        // remove sigular edges
        if ((node->status & REMOVED) == 0 &&
            (node->marked & MARKED) != 0)
        {
            for (k = 0; k < node->max_F; k++)
            {
                edge = &(node->edges[k + node->max_B]);
                d_node = &(graph->nodes[edge->dest]);               
                if (edge->type != INVALID &&
                    (d_node->status & REMOVED) != 0)
                {                  
                    edge->type = INVALID;
                    node->num_F--;
                    node->status |= CHANGED;
                }
            }
            for (k = 0; k < node->max_B; k++)
            {
                edge = &(node->edges[k]);
                d_node = &(graph->nodes[edge->dest]);
                if (edge->type != INVALID &&
                    (d_node->status & REMOVED) != 0)
                {                  
                    edge->type = INVALID;
                    node->num_B--;
                    node->status |= CHANGED;
                }
            }
            node->marked &= SUSPECT;
        }
    }
     
    // reset status
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {
        node_t *node;
        node = &(graph->nodes[i]);
        char old_B_type;
        char old_F_type;
        // reset status
        if ((node->status & REMOVED) == 0 &&
            (node->status & CHANGED) != 0)
        {
            old_B_type = node->B_type;
            old_F_type = node->F_type;
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
                node->status = REMOVED;
                continue;
            }
            // remark
            node->status &= REMOVED;

            if (old_F_type != node->F_type || old_B_type != node->B_type)
            {
                node->status |= CHANGED;  
            }
        }
    }
#ifdef DEBUG2
    fclose(fp);
#endif
}


static void
init_path (string_graph_t * graph)
{
    INDEX_TYPE_T i;    
    
    /* extend path */
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {
        int k;
        node_t *node;
        edge_t *edge;
        path_q_t *path_q;
        path_t *path;
        node = &(graph->nodes[i]);
        path_q = node->path_q;
        if ((node->status & REMOVED) == 0 &&
            ((node->B_type == BRANCH && node->F_type != ENDING) ||
                 (node->F_type == BRANCH && node->B_type != ENDING)))
        {
            for (k = 0; k < node->max_B; k++)
            {
                edge = &(node->edges[k]);
                path = &(path_q->path[k]);
                if (edge->type != INVALID)
                {
                    extend_path (path, graph);
                }
            }
            for (k = 0; k < node->max_F; k++)
            {
                edge = &(node->edges[k + node->max_B]);
                path = &(path_q->path[k + node->max_B]);
                if (edge->type != INVALID)
                {
                    extend_path (path, graph);
                }
            }
        }
    }
}


static void
remove_tips (string_graph_t * graph, int max_tiplen, int num_iters)
{
    INDEX_TYPE_T i;
    int finished;
    int read_len;
    int tiplen;
    int stride;

    finished = 0;
    read_len = graph->read_len;
    tiplen = 1;
    stride = max_tiplen / num_iters;

    while (finished == 0 || num_iters > 0)
    {
        if (finished == 1)
        {
            tiplen += stride;
            num_iters--;
            tiplen = (tiplen > max_tiplen ? max_tiplen : tiplen);
        }
        else
        {
            finished = 1;
        }

        /* extend path */
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            int k;
            node_t *node;
            node_t *d_node;
            edge_t *edge;
            path_q_t *path_q;
            path_t *path;
            node = &(graph->nodes[i]);
            path_q = node->path_q;

            if ((node->status & REMOVED) == 0 &&
                ((node->B_type == BRANCH && node->F_type != ENDING) ||
                 (node->F_type == BRANCH && node->B_type != ENDING)))
            {
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    d_node = &(graph->nodes[path->end]);
                    if (edge->type != INVALID)
                    {
                        if ((d_node->status & CHANGED2) != 0)
                        {
                            path->end = edge->dest;
                            path->hop = edge->dest;
                            path->len = 1;
                            path->end_len = graph->read_len - edge->overlap_len;
                            path->bubble_len = 0;
                            path->type = NONE;
                            path->direction = (edge->type == BB ? 1 : 0);
                            extend_path (path, graph);
                        }
                        else if ((d_node->status & CHANGED) != 0)
                        {
                            extend_path (path, graph);
                        }
                    }
                }
                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    path = &(path_q->path[k + node->max_B]);
                    d_node = &(graph->nodes[path->end]);
                    if (edge->type != INVALID)
                    {
                        if ((d_node->status & CHANGED2) != 0)
                        {
                            path->end = edge->dest;
                            path->hop = edge->dest;
                            path->len = 1;
                            path->end_len = graph->read_len - edge->overlap_len;
                            path->bubble_len = 0;
                            path->type = NONE;
                            path->direction = (edge->type == FB ? 1 : 0);
                            extend_path (path, graph);
                        }
                        else if ((d_node->status & CHANGED) != 0)
                        {
                            extend_path (path, graph);
                        }
                    }
                }
            }
        }

        /* clear flags */
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            node = &(graph->nodes[i]);
            node->status &= REMOVED;
        }

        /* Find tips */
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            int k;
            node_t *node;
            node_t *e_node;
            edge_t *edge;
            path_q_t *path_q;
            path_t *path;            
            int direction;
            int end_len;
            int len;
            node = &(graph->nodes[i]);
            path_q = node->path_q;
            // from a branch node
            if ((node->status & REMOVED) == 0 &&
                ((node->B_type == BRANCH && node->F_type != ENDING) ||
                 (node->F_type == BRANCH && node->B_type != ENDING)))
            {
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    if (edge->type == INVALID ||
                        path->updated != 1 || path->type != TIP)
                    {
                        continue;
                    }
                    e_node = &(graph->nodes[path->end]);
                    len = e_node->path_q == NULL ?
                        0 : (path->direction == 0 ?
                             e_node->path_q->max_B_len :
                             e_node->path_q->max_F_len);
                    len += path->len;
                    end_len = e_node->path_q == NULL ?
                        0 : (path->direction == 0 ?
                             e_node->path_q->max_B_endlen :
                             e_node->path_q->max_F_endlen);
                    end_len += path->end_len - (read_len - edge->overlap_len);
                    if (end_len < max_tiplen && len <= tiplen)
                    {
#ifdef DEBUG2
                        direction = (edge->type == BB ? 1 : 0);
                        remove_path3 (edge->dest, direction, graph);
#else
                        direction = path->direction;
                        remove_path4 (edge->dest, path, direction, graph);
#endif
                        edge->type = INVALID;
                        node->num_B--;
                        len = (len > 63535 ? 63535 : len);
                        end_len = (end_len > 63535 ? 63535 : end_len);
                        if (len > path_q->max_B_len)
                        {
                            path_q->max_B_len = len;
                        }
                        if ((end_len + read_len - edge->overlap_len) >
                            path_q->max_B_endlen)
                        {
                            path_q->max_B_endlen =
                                end_len + read_len - edge->overlap_len;
                        }
                        node->status |= CHANGED;
                    }
                    else if (end_len >= read_len || len > max_tiplen)
                    {
                        path->updated = 0;
                    }
                }
                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    path = &(path_q->path[k + node->max_B]);
                    if (edge->type == INVALID ||
                        path->updated != 1 || path->type != TIP)
                    {
                        continue;
                    }
                    e_node = &(graph->nodes[path->end]);
                    len = e_node->path_q == NULL ?
                        0 : (path->direction == 0 ?
                             e_node->path_q->max_B_len :
                             e_node->path_q->max_F_len);
                    len += path->len;
                    end_len = e_node->path_q == NULL ?
                        0 : (path->direction == 0 ?
                             e_node->path_q->max_B_endlen :
                             e_node->path_q->max_F_endlen);
                    if (end_len < max_tiplen && len <= tiplen)
                    {
#ifdef DEBUG2
                        direction = (edge->type == FB ? 1 : 0);
                        remove_path3 (edge->dest, direction, graph);
#else
                        direction = path->direction;
                        remove_path4 (edge->dest, path, direction, graph);
#endif
                        edge->type = INVALID;
                        node->num_F--;
                        end_len += read_len - edge->overlap_len;
                        len = (len > 63535 ? 63535 : len);
                        end_len = (end_len > 63535 ? 63535 : end_len);
                        if (len > path_q->max_F_len)
                        {
                            path_q->max_F_len = len;
                        }
                        if (end_len > path_q->max_F_endlen)
                        {
                            path_q->max_F_endlen = end_len;
                        }
                        node->status |= CHANGED;
                    }
                    else if (end_len >= read_len || len > max_tiplen)
                    {
                        path->updated = 0;
                    }
                }
            }
        }

        /* set status */
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            node = &(graph->nodes[i]);
            char old_B_type;
            char old_F_type;
            // reset status
            if ((node->status & REMOVED) == 0 &&
                (node->status & CHANGED) != 0)
            {
                old_B_type = node->B_type;
                old_F_type = node->F_type;
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
                    node->status = REMOVED;
                    continue;
                }
                // remark
                node->status &= REMOVED;

                if (old_F_type != node->F_type || old_B_type != node->B_type)
                {
                    if ((node->F_type != BRANCH && node->B_type != BRANCH) ||
                        node->F_type == ENDING || node->B_type == ENDING)
                    {
                        if (node->path_q->path != NULL)
                        {
                            free (node->path_q->path);
                            node->path_q->path = NULL;
                        }
                    }
                    node->status |= CHANGED;
                    finished = 0;
                }
            }
        }
    }

#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {
        int k;
        node_t *node;
        edge_t *edge;
        node_t *d_node;
        path_q_t * e_path_q;
        int end_len;
        
        node = &(graph->nodes[i]); 
        // remove sigular edges
        if ((node->status & REMOVED) == 0 && (node->marked & MARKED) != 0)
        {
            if ((node->B_type == BRANCH && node->F_type != ENDING) ||
                (node->F_type == BRANCH && node->B_type != ENDING))
            {
                continue;
            }

            for (k = 0; k < node->max_F; k++)
            {
                edge = &(node->edges[k + node->max_B]);
                d_node = &(graph->nodes[edge->dest]);
                if (edge->type != INVALID &&
                    (d_node->status & REMOVED) != 0)
                {
                    if (node->path_q == NULL)
                    {
                        node->path_q = (path_q_t *)malloc (sizeof(path_q_t));
                        assert (node->path_q != NULL);
                        node->path_q->max_F_endlen = 0;
                        node->path_q->max_B_endlen = 0;
                    }
                    e_path_q = d_node->path_q;
                    end_len = e_path_q == NULL ?
                              0 : (edge->type == FF ?
                                   e_path_q->max_B_endlen :
                                   e_path_q->max_F_endlen);
                    end_len += read_len - edge->overlap_len;
                    if (end_len > node->path_q->max_F_endlen)
                    {
                        node->path_q->max_F_endlen = end_len;
                    }
                    edge->type = INVALID;
                    node->num_F--;
                }
            }
            for (k = 0; k < node->max_B; k++)
            {
                edge = &(node->edges[k]);
                d_node = &(graph->nodes[edge->dest]);
                if (edge->type != INVALID &&
                    (d_node->status & REMOVED) != 0)
                {
                    if (node->path_q == NULL)
                    {
                        node->path_q = (path_q_t *)malloc (sizeof(path_q_t));
                        assert (node->path_q != NULL);
                        node->path_q->max_F_endlen = 0;
                        node->path_q->max_B_endlen = 0;
                    }
                    e_path_q = d_node->path_q;
                    end_len = e_path_q == NULL ?
                              0 : (edge->type == BF ?
                                   e_path_q->max_B_endlen :
                                   e_path_q->max_F_endlen);
                    end_len += read_len - edge->overlap_len;
                    if (end_len > node->path_q->max_B_endlen)
                    {
                        node->path_q->max_B_endlen = end_len;
                    }
                    edge->type = INVALID;
                    node->num_B--;
                }
            }
            
            node->marked &= SUSPECT;
            
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
                node->status = REMOVED;
            }
        }       
    }
}


// remove bubbles
static void
remove_bubbles (string_graph_t * graph, int min_bubblelen)
{
#ifdef DEBUG2
    FILE *fp;
    fp = fopen ("bubble.dat", "w+");
#endif

    INDEX_TYPE_T i;
    int read_len;
    int min_overlap;
    double average_overlap;
    int finished;

    read_len = graph->read_len;
    min_overlap = graph->min_overlap;
    average_overlap = graph->average_overlap;
    finished = 0;
    
    while (finished == 0)
    {
        finished = 1;

        // extend path
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            int k;
            node_t *node;
            node_t *d_node;
            edge_t *edge;
            path_q_t *path_q;
            path_t *path;
            node = &(graph->nodes[i]);
            path_q = node->path_q;
        
            if ((node->status & REMOVED) == 0 &&
                ((node->B_type == BRANCH && node->F_type != ENDING) ||
                 (node->F_type == BRANCH && node->B_type != ENDING)))
            {
                path_q->max_B_len = 0;
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    d_node = &(graph->nodes[path->end]);
                    if (edge->type != INVALID)
                    {
                        if ((d_node->status & CHANGED2) != 0)
                        {
                            path->end = edge->dest;
                            path->hop = edge->dest;
                            path->len = 1;
                            path->end_len = graph->read_len - edge->overlap_len;
                            path->bubble_len = 0;
                            path->type = NONE;
                            path->direction = (edge->type == BB ? 1 : 0);
                            extend_path (path, graph);
                        }
                        else if ((d_node->status & CHANGED) != 0)
                        {
                            extend_path (path, graph);
                        }
                        if (edge->overlap_len > path_q->max_B_len)
                        {
                            path_q->max_B_len = edge->overlap_len;
                        }
                    }
                }
                path_q->max_F_len = 0;
                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    path = &(path_q->path[k + node->max_B]);
                    d_node = &(graph->nodes[path->end]);                   
                    if (edge->type != INVALID)
                    {
                        if ((d_node->status & CHANGED2) != 0)
                        {
                            path->end = edge->dest;
                            path->hop = edge->dest;
                            path->len = 1;
                            path->end_len = graph->read_len - edge->overlap_len;
                            path->bubble_len = 0;
                            path->type = NONE;
                            path->direction = (edge->type == FB ? 1 : 0);
                            extend_path (path, graph);
                        }
                        else if ((d_node->status & CHANGED) != 0)
                        {
                            extend_path (path, graph);
                        }                    
                        if (edge->overlap_len > path_q->max_F_len)
                        {
                            path_q->max_F_len = edge->overlap_len;
                        }
                    }
                }
            }
        }
        
        /* clear flags */
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            node = &(graph->nodes[i]);
            node->status &= REMOVED;        
        }

        // remove bubbles
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            node_t *e_node;
            edge_t *edge;
            edge_t *n_edge;
            path_q_t *path_q;
            path_q_t *e_path_q;
            path_t *path;
            path_t *n_path;
            int k;
            int j;
            int match_len;
            int max_overlaps;
            int bubble_len;
            int e_max_overlaps;
            double depth;
#ifdef DEBUG2
            int direction;
#else
            node_t *n_node;
#endif

            node = &(graph->nodes[i]);
            path_q = node->path_q;
            /* B end */
            if ((node->status & REMOVED) == 0 &&
                node->B_type == BRANCH &&
                node->F_type != ENDING)
            {
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    if (edge->type == INVALID || path->type != BUBBLE)
                    {
                        continue;
                    }
                    e_node = &(graph->nodes[path->end]);
                    e_path_q = e_node->path_q;
                    max_overlaps = read_len - path_q->max_B_len - 1;
                    e_max_overlaps = read_len - 1 -
                        (path->direction == 1 ?
                         e_path_q->max_B_len : e_path_q->max_F_len);
                    match_len = read_len - path->end_len +
                        max_overlaps + e_max_overlaps;
                    depth =
                        (double) read_len -
                        (double) path->end_len / path->len;
                    bubble_len =
                        path->bubble_len - (read_len - edge->overlap_len);
                    if (bubble_len < read_len && match_len < min_overlap)
                    {
                        if (depth < average_overlap)
                        {
#ifdef DEBUG2
                            fprintf (fp, "%d\n", edge->dest);
#endif

#ifdef DEBUG2
                            direction = (edge->type == BB ? 1 : 0);
                            remove_path0 (edge->dest, direction, graph);
#else                          
                            n_node = &(graph->nodes[edge->dest]);
                            if ((edge->type == BF && n_node->F_type == ONE) ||
                                (edge->type == BB && n_node->B_type == ONE))
                            {
                                n_node->status = REMOVED;
                            }
#endif                            
                            edge->type = INVALID;
                            node->num_B--;
                            node->status |= CHANGED;
                        }
                        else
                        {
                            for (j = 0; j < node->max_B; j++)
                            {
                                n_edge = &(node->edges[j]);
                                n_path = &(path_q->path[j]);
                                if (n_edge->type != INVALID &&
                                    j != k &&
                                    n_path->type == BUBBLE &&
                                    n_path->direction == path->direction &&
                                    n_path->end == path->end &&
                                    path->len + min_bubblelen < n_path->len)
                                {
#ifdef DEBUG2
                                   fprintf (fp, "%d\n", edge->dest);
#endif                                
#ifdef DEBUG2
                                   direction = (edge->type == BB ? 1 : 0);
                                   remove_path0 (edge->dest, direction, graph);
#else                                
                                    n_node = &(graph->nodes[edge->dest]);
                                    if ((edge->type == BF
                                         && n_node->F_type == ONE)
                                        || (edge->type == BB
                                            && n_node->B_type == ONE))
                                    {
                                        n_node->status = REMOVED;
                                    }
#endif                                    
                                    edge->type = INVALID;
                                    node->num_B--;
                                    node->status |= CHANGED;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            /* F end */
            if ((node->status & REMOVED) == 0 &&
                node->F_type == BRANCH &&
                node->B_type != ENDING)
            {
                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    path = &(path_q->path[k + node->max_B]);
                    if (edge->type == INVALID || path->type != BUBBLE)
                    {
                        continue;
                    }
                    e_node = &(graph->nodes[path->end]);
                    e_path_q = e_node->path_q;
                    max_overlaps = read_len - path_q->max_F_len - 1;
                    e_max_overlaps = read_len - 1 -
                        (path->direction == 1 ?
                         e_path_q->max_B_len : e_path_q->max_F_len);
                    match_len = read_len - path->end_len +
                        max_overlaps + e_max_overlaps;
                    depth =
                        (double) read_len -
                        (double) path->end_len / path->len;
                    bubble_len =
                        path->bubble_len - (read_len - edge->overlap_len);
                    if (bubble_len < read_len && match_len < min_overlap)
                    {
                        if (depth < average_overlap)
                        {
#ifdef DEBUG2
                            fprintf (fp, "%d\n", edge->dest);
#endif

#ifdef DEBUG2
                            direction = (edge->type == FB ? 1 : 0);
                            remove_path0 (edge->dest, direction, graph);
#else                         
                            n_node = &(graph->nodes[edge->dest]);
                            if ((edge->type == BF && n_node->F_type == ONE) ||
                                (edge->type == BB && n_node->B_type == ONE))
                            {
                                n_node->status = REMOVED;
                            }
#endif                            
                            edge->type = INVALID;
                            node->num_F--;
                            node->status |= CHANGED;
                        }
                        else
                        {
                            for (j = 0; j < node->max_F; j++)
                            {
                                n_edge = &(node->edges[j + node->max_B]);
                                n_path = &(path_q->path[j + node->max_B]);
                                if (n_edge->type != INVALID &&
                                    j != k &&
                                    n_path->type == BUBBLE &&
                                    n_path->direction == path->direction &&
                                    n_path->end == path->end &&
                                    path->len + min_bubblelen < n_path->len)
                                {
#ifdef DEBUG2
                                   fprintf (fp, "%d\n", edge->dest);
#endif

#ifdef DEBUG2
                                   direction = (edge->type == FB ? 1 : 0);
                                   remove_path0 (edge->dest, direction, graph);
#else
                                    n_node = &(graph->nodes[edge->dest]);
                                    if ((edge->type == FF
                                         && n_node->F_type == ONE)
                                        || (edge->type == FB
                                            && n_node->B_type == ONE))
                                    {
                                        n_node->status = REMOVED;
                                    }
#endif                                    
                                    edge->type = INVALID;
                                    node->num_F--;
                                    node->status |= CHANGED;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // reset status
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            char old_B_type;
            char old_F_type;
        
            node = &(graph->nodes[i]);   
            if ((node->status & CHANGED) != 0)
            {
                old_B_type = node->B_type;
                old_F_type = node->F_type;
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
                    node->status = REMOVED;
                    continue;
                }
                node->status &= REMOVED;

                if (old_F_type != node->F_type || old_B_type != node->B_type)
                {
                    if ((node->F_type != BRANCH && node->B_type != BRANCH) ||
                        node->F_type == ENDING || node->B_type == ENDING)
                    {
                        if (node->path_q->path != NULL)
                        {
                            free (node->path_q->path);
                            node->path_q->path = NULL;
                        }
                    }
                    node->status |= CHANGED;
                    finished = 0;
                }
            }
        }
    }
#ifdef DEBUG2
    fclose (fp);
#endif
}


// remove bubble combo
static void
remove_bubblecombo (string_graph_t * graph, int min_bubblelen)
{
#ifdef DEBUG2
    FILE *fp;
    fp = fopen ("junc.dat", "w+");
#endif

    INDEX_TYPE_T i;
    int read_len;
    int min_overlap;
    double average_overlap;
    int finished;

    read_len = graph->read_len;
    min_overlap = graph->min_overlap;
    average_overlap = graph->average_overlap;    
    finished = 0;
    
    while (finished == 0)
    {
        finished = 1;

        // extend path
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            int k;
            node_t *node;
            node_t *d_node;
            edge_t *edge;
            path_q_t *path_q;
            path_t *path;
            node = &(graph->nodes[i]);
            path_q = node->path_q;

            if ((node->status & REMOVED) == 0 &&
                ((node->B_type == BRANCH && node->F_type != ENDING) ||
                (node->F_type == BRANCH && node->B_type != ENDING)))
            {
                path_q->max_B_len = 0;
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    d_node = &(graph->nodes[path->end]);
                    if (edge->type != INVALID)
                    {
                        if ((d_node->status & CHANGED2) != 0)
                        {
                            path->end = edge->dest;
                            path->hop = edge->dest;
                            path->len = 1;
                            path->end_len = graph->read_len - edge->overlap_len;
                            path->bubble_len = 0;
                            path->type = NONE;
                            path->direction = (edge->type == BB ? 1 : 0);
                            extend_path (path, graph);
                        }
                        else if ((d_node->status & CHANGED) != 0)
                        {
                            extend_path (path, graph);
                        }                   
                        if (edge->overlap_len > path_q->max_B_len)
                        {
                            path_q->max_B_len = edge->overlap_len;
                        }
                    }
                }
                path_q->max_F_len = 0;
                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    path = &(path_q->path[k + node->max_B]);
                    d_node = &(graph->nodes[path->end]);                   
                    if (edge->type != INVALID)
                    {
                        if ((d_node->status & CHANGED2) != 0)
                        {
                            path->end = edge->dest;
                            path->hop = edge->dest;
                            path->len = 1;
                            path->end_len = graph->read_len - edge->overlap_len;
                            path->bubble_len = 0;
                            path->type = NONE;
                            path->direction = (edge->type == FB ? 1 : 0);
                            extend_path (path, graph);
                        }
                        else if ((d_node->status & CHANGED) != 0)
                        {
                            extend_path (path, graph);
                        }                    
                        if (edge->overlap_len > path_q->max_F_len)
                        {
                            path_q->max_F_len = edge->overlap_len;
                        }
                    }
                }
            }
        }

        /* clear flags */
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            node = &(graph->nodes[i]);
            node->status &= REMOVED;
        }
        
        // remove junctions
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            int k;
            int j;
            int m;
            path_t *e_path;
            edge_t *e_edge;
            edge_t *edge;
            path_q_t *path_q;
            int flag;
            path_t *path;
            path_t *path1;
            path_t *path2;
            edge_t *edge1;
            edge_t *edge2;
            node_t *e_node1;
            node_t *e_node2;
            path_q_t *e_path_q1;
            path_q_t *e_path_q2;
            int bubble_len1;
            int bubble_len2;
            int end_len1;
            int end_len2;
            int len1;
            int len2;
            int len;
            int e_max_overlaps1;
            int e_max_overlaps2;
            double depth;
            int match_len;
            int bubblelen;
            int num_F;
            int num_B;
#ifdef DEBUG2
            int direction;
#endif

            node = &(graph->nodes[i]);
            path_q = node->path_q;

            if ((node->status & REMOVED) == 0 &&
                ((node->B_type == BRANCH && node->F_type != ENDING) ||
                (node->F_type == BRANCH && node->B_type != ENDING)))
            {
                flag = 0;
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    path = &(path_q->path[k]);
                    path->updated = 0;
                    if (edge->type != INVALID && path->type != BUBBLE)
                    {
                        flag = 1;
                        break;
                    }
                }

                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    path = &(path_q->path[k + node->max_B]);
                    path->updated = 0;
                    if (edge->type != INVALID && path->type != BUBBLE)
                    {
                        flag = 1;
                        break;
                    }
                }
                // junction
                if (flag != 1)
                {
                    for (k = 0; k < node->max_B; k++)
                    {
                        edge1 = &(node->edges[k]);
                        path1 = &(path_q->path[k]);
                        if (edge1->type == INVALID)
                        {
                            continue;
                        }
                        e_node1 = &(graph->nodes[path1->end]);
                        e_path_q1 = e_node1->path_q;
                        e_max_overlaps1 = read_len - 1 -
                            (path1->direction == 1 ?
                             e_path_q1->max_B_len : e_path_q1->max_F_len);
                        bubble_len1 = path1->bubble_len;
                        len1 = path1->len;
                        end_len1 = path1->end_len;
                        for (j = 0; j < node->max_F; j++)
                        {
                            edge2 = &(node->edges[j + node->max_B]);
                            path2 = &(path_q->path[j + node->max_B]);
                            if (edge2->type == INVALID)
                            {
                                continue;
                            }
                            e_node2 = &(graph->nodes[path2->end]);
                            e_path_q2 = e_node2->path_q;
                            e_max_overlaps2 = read_len - 1 -
                                (path2->direction == 1 ?
                                 e_path_q2->max_B_len : e_path_q2->max_F_len);
                            bubble_len2 = path2->bubble_len;
                            len2 = path2->len;
                            end_len2 = path2->end_len;
                            // examine a pair
                            bubblelen = bubble_len1 + bubble_len2;
                            match_len = read_len - (end_len1 + end_len2) +
                                e_max_overlaps1 + e_max_overlaps2;
                            depth = (double) read_len -
                                (double) (end_len1 + end_len2) / (len1 + len2);
                            
                            if (bubblelen < read_len && match_len < min_overlap)
                            {
                                if (depth < average_overlap)
                                {
                                    path1->updated++;
                                    path2->updated++;
                                }
                                else
                                {
                                    len = len1 + len2;
                                    if (path2->direction == 1)
                                    {
                                        for (m = 0; m < e_node2->max_B; m++)
                                        {
                                            e_edge = &(e_node2->edges[m]);
                                            e_path = &(e_path_q2->path[m]);
                                            if (e_edge->type != INVALID &&
                                                e_edge->dest != path2->hop &&
                                                e_edge->dest != i &&
                                                e_path->type == BUBBLE &&
                                                e_path->direction == path1->direction &&
                                                e_path->end == path1->end &&
                                                len + min_bubblelen < e_path->len)
                                            {                                           
                                                path1->updated++;
                                                path2->updated++;
                                                break;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        for (m = 0; m < e_node2->max_F; m++)
                                        {
                                            e_edge = &(e_node2->edges[m + e_node2->max_B]);
                                            e_path = &(e_path_q2->path[m + e_node2->max_B]);
                                            if (e_edge->type != INVALID &&
                                                e_edge->dest != path2->hop &&
                                                e_edge->dest != i &&
                                                e_path->type == BUBBLE &&
                                                e_path->direction == path1->direction &&
                                                e_path->end == path1->end &&
                                                len + min_bubblelen < e_path->len)
                                            {                                           
                                                path1->updated++;
                                                path2->updated++;
                                                break;
                                            }
                                        }
                                        
                                    }
                                }
                            }
                        }
                    }
                    num_F = node->num_F;
                    num_B = node->num_B;
                    for (k = 0; k < node->max_B; k++)
                    {
                        edge = &(node->edges[k]);
                        path = &(path_q->path[k]);
                        if (edge->type != INVALID &&
                            path->updated == num_F)
                        {
#ifdef DEBUG2
                            direction = (edge->type == BB ? 1 : 0);
                            remove_path1 (edge->dest, direction, graph);
#else
                            remove_path2 (edge->dest, path, graph);
#endif
                            node->status |= CHANGED;
                            edge->type = INVALID;
                            node->num_B--;
                        }
                    }
                    for (k = 0; k < node->max_F; k++)
                    {
                        edge = &(node->edges[k + node->max_B]);
                        path = &(path_q->path[k + node->max_B]);
                        if (edge->type != INVALID &&
                            path->updated == num_B)
                        {
#ifdef DEBUG2
                            direction = (edge->type == FB ? 1 : 0);
                            remove_path1 (edge->dest, direction, graph);
#else
                            remove_path2 (edge->dest, path, graph);
#endif
                            node->status |= CHANGED;
                            edge->type = INVALID;
                            node->num_F--;
                        }
                    }
#ifdef DEBUG2       
                    if ((node->status & CHANGED) != 0)
                        fprintf (fp, "%d\n", i);
#endif
                }
            }
        }

        // reset status
#ifdef __OPENMP__
        #pragma omp parallel for
#endif
        for (i = 0; i < graph->num_reads; i++)
        {
            node_t *node;
            node_t *n_node;
            int k;
            int j;
            edge_t *edge;
            edge_t *n_edge;
            char old_B_type;
            char old_F_type;
            int flag;

            // reset status
            node = &(graph->nodes[i]);
            if ((node->status & REMOVED) == 0 && (node->marked & MARKED) != 0)
            {
                for (k = 0; k < node->max_F; k++)
                {
                    edge = &(node->edges[k + node->max_B]);
                    if (edge->type == INVALID)
                    {
                        continue;
                    }
                    n_node = &(graph->nodes[edge->dest]);
                    if ((n_node->status & REMOVED) != 0)
                    {
                        edge->type = INVALID;
                        node->num_F--;
                        node->status |= CHANGED;
                    }
                    else if ((n_node->status & CHANGED) != 0)
                    {
                        flag = 0;
                        // has a FF edge
                        if (edge->type == FF)
                        {
                            // check whether n_node has a FF edge
                            for (j = 0; j < n_node->max_F; j++)
                            {
                                n_edge = &(n_node->edges[j + n_node->max_B]);
                                if (n_edge->type == FF &&
                                    n_edge->dest == i &&
                                    n_edge->overlap_len == edge->overlap_len)
                                {
                                    flag = 1;
                                    break;
                                }
                            }

                        }
                        // has a FB edge
                        else
                        {
                            // check whether n_node has a BF edge
                            for (j = 0; j < n_node->max_B; j++)
                            {
                                n_edge = &(n_node->edges[j]);
                                if (n_edge->type == BF &&
                                    n_edge->dest == i &&
                                    n_edge->overlap_len == edge->overlap_len)
                                {
                                    flag = 1;
                                    break;
                                }
                            }
                        }
                        if (flag == 0)
                        {
                            edge->type = INVALID;
                            node->num_F--;
                            node->status |= CHANGED;
                        }
                    }
                }
                for (k = 0; k < node->max_B; k++)
                {
                    edge = &(node->edges[k]);
                    if (edge->type == INVALID)
                    {
                        continue;
                    }
                    n_node = &(graph->nodes[edge->dest]);
                    if ((n_node->status & REMOVED) != 0)
                    {
                        edge->type = INVALID;
                        node->num_B--;
                        node->status |= CHANGED;
                    }
                    else if ((n_node->status & CHANGED) != 0)
                    {
                        flag = 0;
                        // has a BF edge
                        if (edge->type == BF)
                        {
                            // check whether n_node has a FB edge
                            for (j = 0; j < n_node->max_F; j++)
                            {
                                n_edge = &(n_node->edges[j + n_node->max_B]);
                                if (n_edge->type == FB &&
                                    n_edge->dest == i &&
                                    n_edge->overlap_len == edge->overlap_len)
                                {
                                    flag = 1;
                                    break;
                                }
                            }

                        }
                        // has a BB edge
                        else
                        {
                            // check whether n_node has a BB edge
                            for (j = 0; j < n_node->max_B; j++)
                            {
                                n_edge = &(n_node->edges[j]);
                                if (n_edge->type == BB &&
                                    n_edge->dest == i &&
                                    n_edge->overlap_len == edge->overlap_len)
                                {
                                    flag = 1;
                                    break;
                                }
                            }
                        }
                        if (flag == 0)
                        {
                            edge->type = INVALID;
                            node->num_B--;
                            node->status |= CHANGED;
                        }
                    }
                }
                node->marked &= SUSPECT;
            }

            if ((node->status & CHANGED) != 0)
            {
                old_B_type = node->B_type;
                old_F_type = node->F_type;
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
                    node->status = REMOVED;
                }
                
                node->status &= REMOVED;

                if (old_F_type != node->F_type || old_B_type != node->B_type)
                {
                    if ((node->F_type != BRANCH && node->B_type != BRANCH) ||
                        node->F_type == ENDING || node->B_type == ENDING)
                    {
                        if (node->path_q->path != NULL)
                        {
                            free (node->path_q->path);
                            node->path_q->path = NULL;
                        }
                    }
                    node->status |= CHANGED;
                    finished = 0;
                }
            }           
        }
    }
#ifdef DEBUG2
    fclose(fp);
#endif
}


// clean short chains
static void
clean_shortchain (string_graph_t * graph, int max_tiplen)
{
    INDEX_TYPE_T i;
    int read_len = graph->read_len;

    /* dead ends */
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {
        node_t *node;
        node_t *n_node;
        edge_t *edge;
        int o_len;
        int k;
        int j;
        int direction;
        int end_len;
        int n_len;
        node = &(graph->nodes[i]);
        if ((node->status & REMOVED) != 0)
        {
            continue;
        }
        
        /* if node B is ENDING */
        if (node->B_type == ENDING)
        {
            if (node->path_q != NULL)
            {
                n_len = node->path_q->max_B_endlen;
            }
            else
            {
                node->path_q = (path_q_t *) malloc (sizeof (path_q_t));
                assert (node->path_q != NULL);
                n_len = 0;
            }
            for (k = 0; k < node->max_F; k++)
            {
                edge = &(node->edges[k + node->max_B]);
                if (edge->type != INVALID)
                {
                    end_len = n_len;
                    o_len = end_len + read_len - edge->overlap_len;
                    n_node = &(graph->nodes[edge->dest]);
                    
                    if (o_len <= max_tiplen)
                    {
                        end_len = o_len;
                        node->status = REMOVED;
                    }
                    else
                    {
                        if (node->F_type == BRANCH)
                        {
                            n_node->marked |= MARKED;    
                        }
                        continue;
                    }
                    direction = (edge->type == FB ? 1 : 0);
                    while ((n_node->status & REMOVED) == 0 &&
                           n_node->B_type == ONE &&
                           n_node->F_type == ONE)
                    {
                        if (edge->type == BB || edge->type == FB)
                            direction = 1;
                        else
                            direction = 0;
                        
                        if (direction == 1)
                        {
                            for (j = 0; j < n_node->max_F; j++)
                            {
                                edge = &(n_node->edges[j + n_node->max_B]);
                                if (edge->type != INVALID)
                                {
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (j = 0; j < n_node->max_B; j++)
                            {
                                edge = &(n_node->edges[j]);
                                if (edge->type != INVALID)
                                {
                                    break;
                                }
                            }
                        }
                        o_len = end_len + read_len - edge->overlap_len;
                        if (o_len <= max_tiplen)
                        {
                            end_len = o_len;
                            n_node->status = REMOVED;
                        }
                        else
                        {
                            if (n_node->path_q == NULL)
                            {
                                n_node->path_q =
                                    (path_q_t *) malloc (sizeof (path_q_t));
                                assert (n_node->path_q != NULL);
                            }

                            if (direction == 1)
                            {
                                n_node->path_q->max_B_endlen =
                                    read_len - end_len;
                            }
                            else
                            {
                                n_node->path_q->max_F_endlen =
                                    read_len - end_len;
                            }
                            break;
                        }
                        n_node = &(graph->nodes[edge->dest]);                      
                    } /* end while */
                    n_node->marked |= MARKED;
                }
            }

            if ((node->status & REMOVED) == 0)
            {
                node->B_type = TENDING;
                node->path_q->max_B_endlen = read_len - n_len;
            }
            
        }
        /* if node F is ENDING */
        if (node->F_type == ENDING)
        {
            if (node->path_q != NULL)
            {
                n_len = node->path_q->max_F_endlen;
            }
            else
            {
                node->path_q = (path_q_t *) malloc (sizeof (path_q_t));
                assert (node->path_q != NULL);
                n_len = 0;
            }
            for (k = 0; k < node->max_B; k++)
            {
                edge = &(node->edges[k]);
                if (edge->type != INVALID)
                {
                    end_len = n_len;
                    o_len = end_len + read_len - edge->overlap_len;
                    n_node = &(graph->nodes[edge->dest]);
                    if (o_len <= max_tiplen)
                    {
                        end_len = o_len;
                        node->status = REMOVED;
                    }
                    else
                    {
                        if (node->B_type == BRANCH)
                        {
                            n_node->marked |= MARKED;    
                        }
                        continue;
                    }
                    
                    while ((n_node->status & REMOVED) == 0 &&
                           n_node->B_type == ONE &&
                           n_node->F_type == ONE)
                    {
                        if (edge->type == BB || edge->type == FB)
                            direction = 1;
                        else
                            direction = 0;
                        if (direction == 1)
                        {
                            for (j = 0; j < n_node->max_F; j++)
                            {
                                edge = &(n_node->edges[j + n_node->max_B]);
                                if (edge->type != INVALID)
                                {
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (j = 0; j < n_node->max_B; j++)
                            {
                                edge = &(n_node->edges[j]);
                                if (edge->type != INVALID)
                                {
                                    break;
                                }
                            }
                        }
                        o_len = end_len + read_len - edge->overlap_len;
                        if (o_len <= max_tiplen)
                        {
                            end_len = o_len;
                            n_node->status = REMOVED;
                        }
                        else
                        {
                            if (n_node->path_q == NULL)
                            {
                                n_node->path_q =
                                    (path_q_t *) malloc (sizeof (path_q_t));
                                assert (n_node->path_q != NULL);
                            }

                            if (direction == 1)
                            {
                                n_node->path_q->max_B_endlen =
                                    read_len - end_len;
                            }
                            else
                            {
                                n_node->path_q->max_F_endlen =
                                    read_len - end_len;
                            }
                            break;
                        }
                        n_node = &(graph->nodes[edge->dest]);
                    } /* end while */
                    n_node->marked |= MARKED;
                }
            }

            if ((node->status & REMOVED) == 0)
            {
                node->F_type = TENDING;
                node->path_q->max_F_endlen = read_len - n_len;
            }
            
        }
    }

#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {
        node_t *node;
        node_t *n_node;
        edge_t *edge;
        path_q_t * e_path_q;
        int end_len;
        int k;
        node = &(graph->nodes[i]);
        if ((node->status & REMOVED) == 0 && (node->marked & MARKED) != 0)
        {
            for (k = 0; k < node->max_F; k++)
            {
                edge = &(node->edges[k + node->max_B]);
                n_node = &(graph->nodes[edge->dest]);
                if (edge->type != INVALID &&
                    (n_node->status & REMOVED) != 0)
                {
                    if (node->path_q == NULL)
                    {
                        node->path_q = (path_q_t *)malloc (sizeof(path_q_t));
                        assert (node->path_q != NULL);
                        node->path_q->max_F_endlen = 0;
                        node->path_q->max_B_endlen = 0;
                    }
                    e_path_q = n_node->path_q;
                    end_len = e_path_q == NULL ?
                              0 : (edge->type == FF ?
                                   e_path_q->max_B_endlen :
                                   e_path_q->max_F_endlen);
                    end_len += read_len - edge->overlap_len;
                    if (end_len > node->path_q->max_F_endlen)
                    {
                        node->path_q->max_F_endlen = end_len;
                    }
                    edge->type = INVALID;
                    node->num_F--;                    
                }
            }
            for (k = 0; k < node->max_B; k++)
            {
                edge = &(node->edges[k]);                
                n_node = &(graph->nodes[edge->dest]);
                if (edge->type != INVALID &&
                    (n_node->status & REMOVED) != 0)
                {
                    if (node->path_q == NULL)
                    {
                        node->path_q = (path_q_t *)malloc (sizeof(path_q_t));
                        assert (node->path_q != NULL);
                        node->path_q->max_F_endlen = 0;
                        node->path_q->max_B_endlen = 0;
                    }
                    e_path_q = n_node->path_q;
                    end_len = e_path_q == NULL ?
                              0 : (edge->type == BF ?
                                   e_path_q->max_B_endlen :
                                   e_path_q->max_F_endlen);
                    end_len += read_len - edge->overlap_len;
                    if (end_len > node->path_q->max_B_endlen)
                    {
                        node->path_q->max_B_endlen = end_len;
                    }
                    edge->type = INVALID;
                    node->num_B--;                    
                }
            }
            node->marked &= SUSPECT;          
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
                    node->B_type = TENDING;
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
                    node->F_type = TENDING;
                }
            }
            else
            {
                node->status = REMOVED;
            }
        }        
    }
}


static void
list_contigs (string_graph_t * graph, char *reads, INDEX_TYPE_T * num_contigs,
              INDEX_TYPE_T * len_contigs, int min_len)
{
    INDEX_TYPE_T i;
    int j;
    int read_len = graph->read_len;
    FILE *fp = fopen ("contigs.fa", "w+");
#ifdef DEBUG2
    FILE *fp1 = fopen ("contigs1.fa", "w+");
    long last_pos1 = 0;
#endif
    long last_pos = 0;
    INDEX_TYPE_T num = 0;
    INDEX_TYPE_T len = 0;

    /* for each branch nodes, remove the multiple edges. */
#ifdef __OPENMP__
    #pragma omp parallel for
#endif
    for (i = 0; i < graph->num_reads; i++)
    {

        int k;
        node_t *node = &(graph->nodes[i]);
        edge_t *edge;
        INDEX_TYPE_T dest;
        char type;
        node_t *d_node;
        char *d_type;

        if ((node->status & REMOVED) != 0)
        {
            continue;
        }

        // B edge type is BRANCE
        if (node->B_type == BRANCH)
        {
            node->B_type = TENDING;
            // for each B edge
            for (k = 0; k < node->max_B; k++)
            {
                edge = &(node->edges[k]);
                if (edge->type == INVALID)
                {
                    continue;
                }
                dest = edge->dest;
                type = edge->type;
                d_node = &(graph->nodes[dest]);
                d_type = (type == BB ? &(d_node->B_type) : &(d_node->F_type));
                if (*d_type == ONE)
                {
                    *d_type = HALF;
                }
                else if (*d_type == BRANCH || *d_type == HALF)
                {
                    node->B_type = HALF;
                }
            }
        }

        // F edge type is BRANCH
        if (node->F_type == BRANCH)
        {
            node->F_type = TENDING;
            // for each F edge
            for (k = 0; k < node->max_F; k++)
            {
                edge = &(node->edges[k + node->max_B]);
                if (edge->type == INVALID)
                {
                    continue;
                }
                dest = edge->dest;
                type = edge->type;
                d_node = &(graph->nodes[dest]);
                d_type = (type == FB ? &(d_node->B_type) : &(d_node->F_type));
                if (*d_type == ONE)
                {
                    *d_type = HALF;
                }
                else if (*d_type == BRANCH || *d_type == HALF)
                {
                    node->F_type = HALF;
                }
            }
        }
    }

    /* get contigs */
    for (i = 0; i < graph->num_reads; i++)
    {
        node_t *node = &(graph->nodes[i]);
        int B_end = 0;
        int F_end = 0;
        int direction = 0;
        int max = 0;
        INDEX_TYPE_T n = 0;
        int flag = 0;
        edge_t *edge;
        char *read;

        if ((node->status & REMOVED) != 0)
        {
            continue;
        }

        /* starting from each HALF or ENDING node */
        if ((node->status & VISITED) == 0
            && (node->B_type == HALF || node->B_type == TENDING))
        {
            direction = 0;
        }
        else if ((node->status & VISITED) == 0
                 && (node->F_type == HALF || node->F_type == TENDING))
        {
            direction = 1;
        }
        else
        {
            continue;
        }

        len = 0;
        n = i;
        // traverse to get one contig
        while (1)
        {
            // start from ENDING or HALF node
            node->status |= VISITED;
            // B edge side
            if (node->B_type == TENDING)
            {
                if (node->path_q != NULL)
                    B_end = node->path_q->max_B_endlen;
                else
                    B_end = 0;
            }
            else if (node->B_type == ONE && direction == 0)
            {
                B_end = 0;
            }
            // get the longest overlap
            else
            {
                max = 0;
                for (j = 0; j < node->max_B; j++)
                {
                    edge = &(node->edges[j]);
                    if (edge->type != INVALID &&
                        max < edge->overlap_len)
                    {
                        max = edge->overlap_len;
                    }
                }
                B_end = max;
            }

            // F edge side
            if (node->F_type == TENDING)
            {
                if (node->path_q != NULL)
                    F_end = node->path_q->max_F_endlen;
                else
                    F_end = 0;
            }
            else if (node->F_type == ONE && direction == 1)
            {
                F_end = 0;
            }
            // get the longest overlap
            else
            {
                max = 0;
                for (j = 0; j < node->max_F; j++)
                {
                    edge = &(node->edges[j + node->max_B]);
                    if (edge->type != INVALID &&
                        max < edge->overlap_len)
                    {
                        max = edge->overlap_len;
                    }
                }
                F_end = max;
            }
            F_end = read_len - F_end;

            // locate the read, start traversal
            read = reads + read_len * n;
            // forward read
            if (direction == 0)
            {
#ifdef DEBUG2
                if (B_end < F_end && flag == 0)
                {
                    last_pos1 = ftell (fp1);
                    fprintf (fp1, ">contig %d\n", num);
                    fprintf (fp1, "%d, ", n);
                }
#endif
                for (j = B_end; j < F_end; j++)
                {
                    // a new contig
                    if (flag == 0)
                    {
                        // store the postion of file
                        last_pos = ftell (fp);
#if defined(__INDEX_U32__)
                        fprintf (fp, ">contig %d\n", num);
#elif defined(__INDEX_U64__)
                        fprintf (fp, ">contig %lld\n", num);
#else
                        fprintf (fp, ">contig %d\n", num);
#endif

                        num++;
                        flag = 1;
                    }
                    fprintf (fp, "%c", nt_table[(int) read[j]]);
                    len++;
                }

                if (node->F_type == TENDING || node->F_type == HALF)
                {
                    break;
                }
                for (j = 0; j < node->max_F; j++)
                {
                    edge = &(node->edges[node->max_B + j]);
                    if (edge->type != INVALID)
                    {
                        break;
                    }
                }
                n = edge->dest;
                node = &(graph->nodes[n]);
                direction = (edge->type == FB ? 0 : 1);
            }
            // reverse read
            else
            {
#ifdef DEBUG2
                if (B_end < F_end && flag == 0)
                {
                    last_pos1 = ftell (fp1);
                    fprintf (fp1, ">contig %d\n", num);
                    fprintf (fp1, "%d, ", n);
                }
#endif
                for (j = F_end - 1; j >= B_end; j--)
                {
                    // a new contig
                    if (flag == 0)
                    {
                        // store the postion of file
                        last_pos = ftell (fp);
#if defined(__INDEX_U32__)
                        fprintf (fp, ">contig %d\n", num);
#elif defined(__INDEX_U64__)
                        fprintf (fp, ">contig %lld\n", num);
#else
                        fprintf (fp, ">contig %d\n", num);
#endif
                        num++;
                        flag = 1;
                    }
                    fprintf (fp, "%c",
                             nt_table[(int) (ntc_table[(int) read[j]])]);
                    len++;
                }

                if (node->B_type == TENDING || node->B_type == HALF)
                {
                    break;
                }
                for (j = 0; j < node->max_B; j++)
                {
                    edge = &(node->edges[j]);
                    if (edge->type != INVALID)
                    {
                        break;
                    }
                }
                n = edge->dest;
                node = &(graph->nodes[n]);
                direction = (edge->type == BB ? 0 : 1);
            }
#ifdef DEBUG2
            if (flag == 1)
                fprintf (fp1, "%d, ", n);
#endif
        }

#ifdef DEBUG2
        if (flag == 1)
            fprintf (fp1, "end\n");
#endif
        // end of a contig, be discarded when the length is shorter than read length
        if (flag == 1)
        {
            if (len > read_len && len >= min_len)
            {
                fprintf (fp, "\n");
                *len_contigs += len;
            }
            else
            {
                // restore the file
                fseek (fp, last_pos, SEEK_SET);
#ifdef DEBUG2
                fseek (fp1, last_pos1, SEEK_SET);
#endif
                num--;
            }
        }

    }

    // end of file
    *num_contigs = num;
#ifdef DEBUG2
    fprintf (fp1, "; This is end of contig file generated by pasqual.\n");
    fprintf (fp1,
             "; Assemble %d reads with %d length, output %d contigs with total length %d. ",
             graph->num_reads, graph->read_len, num, len);
#endif

    fprintf (fp, "; This is end of contig file generated by pasqual.\n");
#if defined(__INDEX_U32__)
    fprintf (fp,
             "; Assemble %d reads with %d length, output %d contigs with total length %d. ",
             graph->num_reads, graph->read_len, num, len);
#elif defined(__INDEX_U64__)
    fprintf (fp,
             "; Assemble %lld reads with %d length, output %lld contigs with total length %lld. ",
             graph->num_reads, graph->read_len, num, len);
#else
    fprintf (fp,
             "; Assemble %d reads with %d length, output %d contigs with total length %d. ",
             graph->num_reads, graph->read_len, num, len);
#endif

    for (i = 0; i < graph->read_len; i++)
    {
        fprintf (fp, "$");
#ifdef DEBUG2
        fprintf (fp1, "$");
#endif
    }
#ifdef DEBUG2
    fprintf (fp1, "\n");
    fclose (fp1);
#endif
    fprintf (fp, "\n");
    fclose (fp);
}


void
assemble (string_graph_t * graph, char *reads, INDEX_TYPE_T * num_contigs,
          INDEX_TYPE_T * len_contigs, int min_len, int max_tiplen,
          int min_bubblelen, int num_iters)
{
#ifdef __TIMING__
    struct timeval tv1, tv2;
    double time_pass;
    gettimeofday (&tv1, NULL);
#endif
    DPRINTF (3, "Extract contigs\n");

    DPRINTF (3, "  preparing overlap graphs ... ");
    
#ifdef DEBUG2
    check_missing_edges (graph);
#endif

    // remove tips and bubbles
    init_path (graph);

    DPRINTF (3, "done\n");

    DPRINTF (3, "  removing tips ... ");

    remove_tips (graph, max_tiplen, num_iters);
    
    DPRINTF (3, "done\n");

    DPRINTF (3, "  removing bubbles ... ");

    remove_bubbles (graph, min_bubblelen);
    
    remove_bridge (graph);

    remove_bubbles (graph, min_bubblelen);
    
    remove_bubblecombo (graph, min_bubblelen);    
    
    DPRINTF (3, "done\n");

    DPRINTF (3, "  removing short chains ... ");

    clean_shortchain (graph, max_tiplen);

    DPRINTF (3, "done\n");

#ifdef DEBUG2
    FILE *fp;
    INDEX_TYPE_T i;
    fp = fopen ("remain.dat", "w+");
    for (i = 0; i < graph->num_reads; i++)
    {
        node_t *node;
        node = &(graph->nodes[i]);
        if ((node->status & REMOVED) == 0 &&
            (node->B_type != ENDING || node->F_type != ENDING))
        {
            fprintf (fp, "%d\n", i);
        }
    }
    fclose (fp);
#endif

#ifdef DEBUG2
    fp = fopen ("soverlaps.dat", "w+");
    for (i = 0; i < graph->num_reads; i++)
    {
        node_t *node;
        edge_t *edge;
        node_t *n_node;
        int j;
        node = &(graph->nodes[i]);
        fprintf (fp, "read %d: ", i);
        if ((node->status & REMOVED) == 0)
        {
            for (j = 0; j < node->max_B; j++)
            {
                edge = &(node->edges[j]);
                if (edge->type != INVALID)
                {
                    n_node = &(graph->nodes[edge->dest]);
                    if ((n_node->status & REMOVED) == 0)
                    {
                        fprintf (fp, "(%d, ", edge->dest);
                        fprintf (fp, "%d, ", edge->overlap_len);
                        if (edge->type == BF)
                        {
                            fprintf (fp, "BF); ");
                        }
                        else
                        {
                            fprintf (fp, "BB); ");
                        }
                    }
                }
            }
            for (j = 0; j < node->max_F; j++)
            {
                edge = &(node->edges[j + node->max_B]);
                if (edge->type != INVALID)
                {
                    n_node = &(graph->nodes[edge->dest]);
                    if ((node->status & REMOVED) == 0)
                    {
                        fprintf (fp, "(%d, ", edge->dest);
                        fprintf (fp, "%d, ", edge->overlap_len);
                        if (edge->type == FF)
                        {
                            fprintf (fp, "FF); ");
                        }
                        else
                        {
                            fprintf (fp, "FB); ");
                        }
                    }
                }
            }
        }
        else
        {
            fprintf (fp, "REMOVED ");
        }
        fprintf (fp, "end\n");
    }
    fclose (fp);
#endif

    *num_contigs = 0;
    *len_contigs = 0;

    DPRINTF (3, "  extracting contigs ... ");

    // list contigs
    list_contigs (graph, reads, num_contigs, len_contigs, min_len);

    DPRINTF (3, "done\n");

#ifdef __TIMING__
    gettimeofday (&tv2, NULL);
    time_pass =
        (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec -
                                              tv1.tv_usec) / 1000.0;
    DPRINTF (3, "  takes %.3lf ms\n", time_pass);
#endif
}
