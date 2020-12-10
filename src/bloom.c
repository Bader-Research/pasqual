#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bloom.h"
#include "common.h"


#define SETBIT(a, n) (a[n/CHAR_BIT] |= (1<<(n%CHAR_BIT)))
#define GETBIT(a, n) (a[n/CHAR_BIT] & (1<<(n%CHAR_BIT)))


/* -------------------------------------------------------------------------- */
/* definition of hash functions                                               */
/* hash2 ~ hash 4 is not independent                                          */

INDEX_TYPE_T
hash0 (char *key, INDEX_TYPE_T size)
{
    INDEX_TYPE_T h;
    h = 0;
    while (*key)
    {
        h ^= (h << 5) + (h >> 2) + (unsigned char) *key++;
    }
    h = h % size;
    return h;
}


INDEX_TYPE_T
hash1 (char *key, INDEX_TYPE_T size)
{
    INDEX_TYPE_T h;
    h = 0;
    while (*key)
    {
        h = (unsigned char) *key++ + (h << 6) + (h << 16) - h;
    }
    h = h % size;
    return h;
}


INDEX_TYPE_T
hash2 (char *key, INDEX_TYPE_T size)
{
    INDEX_TYPE_T h;
    h = 0;
    while (*key)
    {
        h += *key++;
    }
    h = h % size;
    return h;
}


INDEX_TYPE_T
hash3 (char *key, INDEX_TYPE_T size)
{
    INDEX_TYPE_T h;
    h = 0;
    while (*key)
    {
        h ^= *key++;
    }
    h = h % size;
    return h;
}


INDEX_TYPE_T
hash4 (char *key, INDEX_TYPE_T size)
{
    INDEX_TYPE_T h;
    h = 0;
    while (*key)
    {
        h = 33 * h + *key++;
    }
    h = h % size;
    return h;
}


hashfunc_t func[MAX_FUNC] = { hash0, hash1, hash2, hash3, hash4 };


int
bloom_filter (unsigned char *array, INDEX_TYPE_T size, int nfunc, char *str,
              char *r_str)
{
    unsigned char bits[MAX_FUNC];
    INDEX_TYPE_T hash_value[MAX_FUNC];
    unsigned char r_bits;
    INDEX_TYPE_T r_hash_value;
    int flag;
    int i;
    flag = 0;

    // match forward reads
    for (i = 0; i < nfunc; i++)
    {
        hash_value[i] = func[i] (str, size);
        bits[i] = GETBIT (array, hash_value[i]);
        if (bits[i] == 0)
            flag = 1;
    }
    if (flag == 0)
    {
        goto match;
    }

    // match reverse reads
    for (i = 0; i < nfunc; i++)
    {
        r_hash_value = func[i] (r_str, size);
        r_bits = GETBIT (array, r_hash_value);
        if (r_bits == 0)
        {
            goto set;
        }
    }

  match:
    return 1;
  set:
    for (i = 0; i < nfunc; i++)
    {
        SETBIT (array, hash_value[i]);
    }
    return 0;
}
