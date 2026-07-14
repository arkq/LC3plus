/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef _LC3PLUS_SCRATCH_ALLOCATOR_H
#define _LC3PLUS_SCRATCH_ALLOCATOR_H

#include "defines.h"
#include <stdlib.h>
#include <inttypes.h>

typedef unsigned short lc3_block_size_t;

#define LC3_SCRATCH_MAX_BLOCKS 32

struct lc3_scratch_data
{
    uintptr_t current_block_address;
    lc3_block_size_t current_block_size;
    size_t alignment;
    size_t address_size;
    size_t max_alloc_size;
    size_t current_alloc_size;
    int max_scratch_calculation_only;
#ifdef VERBOSE_SCRATCH_ALLOC
    lc3_block_size_t alloc_sizes[LC3_SCRATCH_MAX_BLOCKS];
    uintptr_t max_address;
    int stack_index;
    int max_stack_index;
#endif
};

typedef struct lc3_scratch_data* lc3_scratch_t;

/* x refers to the number of calls of lc3_scratch_push() */
#ifdef VERBOSE_SCRATCH_ALLOC
#  ifdef _MSC_VER
#    define SCRATCH_BUFFER_OVERHEAD( x ) ( sizeof( struct lc3_scratch_data ) + x * ( 4 + ( ( (size_t) 1 ) << SCRATCH_BUFFER_ALIGNMENT_BITS ) ) )
#  else
#    define SCRATCH_BUFFER_OVERHEAD( x ) ( sizeof( struct lc3_scratch_data ) + x * ( 4 + ( 1 << SCRATCH_BUFFER_ALIGNMENT_BITS ) ) )
#  endif
#else
#  ifdef _MSC_VER
#    define SCRATCH_BUFFER_OVERHEAD( x ) ( sizeof( struct lc3_scratch_data ) + x * ( ( (size_t) 1 ) << SCRATCH_BUFFER_ALIGNMENT_BITS ) )
#  else
#    define SCRATCH_BUFFER_OVERHEAD( x ) ( sizeof( struct lc3_scratch_data ) + x * ( 1 << SCRATCH_BUFFER_ALIGNMENT_BITS ) )
#  endif
#endif

/* int max_scratch_calculation: 1 -> Determine max. used scratch without allocating anything */
lc3_scratch_t lc3_scratch_init( void* start, size_t size, size_t alignment_bits, int max_scratch_calculation );

#if defined( VERBOSE_SCRATCH_ALLOC ) || defined( DEBUG_MEM )
void* _lc3_scratch_push( lc3_scratch_t scratch, size_t size, const char* func, const char* file, int line );
void* _lc3_scratch_pop( lc3_scratch_t scratch, void* ptr, const char* func, const char* file, int line );
#else
void* lc3_scratch_push( lc3_scratch_t scratch, size_t size );
void* lc3_scratch_pop( lc3_scratch_t scratch, void* ptr );
#endif
#ifdef VERBOSE_SCRATCH_ALLOC
void lc3_scratch_print_status( lc3_scratch_t scratch );
#endif

#endif /* SCRATCH_ALLOC */
