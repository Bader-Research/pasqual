#ifndef __COMMON_H__
#define __COMMON_H__

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>


//#define DEBUG1
//#define DEBUG2

/* -------------------------------------------------------------------------- */
/* type of data structures                                                    */

typedef signed long long i64;
typedef unsigned long long u64;
typedef signed int i32;
typedef unsigned int u32;
typedef signed char i8;
typedef unsigned char u8;
typedef unsigned short u16;
typedef signed short i16;

enum Boolean
{ false = 0, true = 1 };
typedef enum Boolean boolean;

// size of index type
#if defined(__INDEX_U32__)

#define INDEX_TYPE_T u32

#elif defined(__INDEX_U64__)

#define  INDEX_TYPE_T u64

#else

#define INDEX_TYPE_T u32

#endif


/* -------------------------------------------------------------------------- */
/* parameters                                                                 */

// number of buckets for intial sort
#define NUM_SLOTS  85
// sort in SLOT_DEPTH runs
#define SLOT_DEPTH    5
// number of buckets per run
#define SLOT_NUM       17
// threshold for insetion sort
#define INSERT_SORT_SIZE 3500
// P, intial parallelism after intial sort
#define SEQ_DEPTH  9


/* -------------------------------------------------------------------------- */
/* debug tools                                                                */

#define _DEBUG_LEVEL_    3      //  0 to 10,
                             //  0 is no debug print info at all,
                             //  10 is full info

#define max(a, b) ((a) > (b) ? (a) : (b))

#define min(a, b) ((a) < (b) ? (a) : (b))

#if ( _DEBUG_LEVEL_ == -1 )
#define DPRINTF( level, fmt, ... )        {}
#else
#define DPRINTF( level, fmt, args... )              \
        do                                                           \
        {                                                             \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ ) \
            {                                                         \
                fprintf( stdout, fmt, ##args );             \
                fflush( stdout );                                 \
            }                                                         \
        } while ( 0 )
#endif

#define ERROR_EXIT(fmt, args... )                    \
    do                                                           \
    {                                                             \
        fprintf (stderr, "Error: ");                          \
        fprintf( stderr, fmt, ##args );                  \
        fflush( stderr );                                      \
        exit (0);                                                \
    } while ( 0 )


/* -------------------------------------------------------------------------- */
/* global variables                                                           */

extern char nt_table[5];
extern char ntc_table[5];


/* -------------------------------------------------------------------------- */
/* definition of functions                                                    */

void *malloc_align (size_t size, unsigned int log2_align);

void free_align (void *ptr);


#endif /* __COMMON_H__ */
