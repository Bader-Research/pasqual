#include "common.h"

#ifdef __OPENMP__
#include <omp.h>
#endif
#ifdef __TIMING__
#include <time.h>
#include <sys/time.h>
#endif


char nt_table[5] = { '@', 'A', 'C', 'G', 'T' };
char ntc_table[5] = { 0, 4, 3, 2, 1 };


// aligned malloc
inline void *
malloc_align (size_t size, unsigned int log2_align)
{
    void *ret;

    unsigned long align = 1 << log2_align;
    char *real = (char *) malloc (size + sizeof (void *) + (align - 1));
    if (real)
    {
        unsigned long offset =
            (align - (unsigned long) (real + sizeof (void *))) & (align - 1);
        ret = (void *) ((real + sizeof (void *)) + offset);
        *((void **) ret - 1) = (void *) real;
    }
    else
    {
        ret = (void *) (real);
    }
    return ret;
}


// aligned realloc
inline void *
realloc_align (void *ptr, size_t size, unsigned int log2_align)
{
    void *ret;
    void *real_ori;

    unsigned long align = 1 << log2_align;
    real_ori = *((void **) (ptr) - 1);
    char *real =
        (char *) realloc (real_ori, size + sizeof (void *) + (align - 1));
    if (real)
    {
        unsigned long offset =
            (align - (unsigned long) (real + sizeof (void *))) & (align - 1);
        ret = (void *) ((real + sizeof (void *)) + offset);
        *((void **) ret - 1) = (void *) real;
    }
    else
    {
        ret = (void *) (real);
    }
    return ret;
}


// aligned free
inline void
free_align (void *ptr)
{
    void *real;

    if (ptr)
    {
        real = *((void **) (ptr) - 1);
        free (real);
    }
}
