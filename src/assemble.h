#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__


#include "string_graph.h"


void assemble (string_graph_t * graph, char *reads,
               INDEX_TYPE_T * num_contigs, INDEX_TYPE_T * len_contigs,
               int min_len, int max_tiplen, int bubblelen, int tip_stride);


#endif /* __ASSEMBLER_H__ */
