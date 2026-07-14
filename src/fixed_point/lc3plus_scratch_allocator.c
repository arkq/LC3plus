/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "lc3plus_scratch_allocator.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

lc3_scratch_t lc3_scratch_init( void* start, size_t size, size_t alignment_bits, int max_scratch_calculation )
{
    size_t alignment;
    size_t offset = 0;
    uintptr_t start_address = (uintptr_t) start;

    /* make sure start address fullfills alignment criteria */
#ifdef _MSC_VER
    alignment = ( (size_t) 1 ) << alignment_bits;
#else
    alignment = 1 << alignment_bits;
#endif
    if ( start_address & ( alignment - 1 ) )
    {
        offset = alignment - ( start_address & ( alignment - 1 ) );
        start_address += offset;

        if ( size < offset )
        {
#ifdef VERBOSE_SCRATCH_ALLOC
            fprintf( stderr, "warning: scratch size %zu too small in lc3_scratch_init", size );
            assert( 0 );
#endif
            return NULL;
        }

        size -= offset;
    }

    /* place lc3_scratch_data struct at beginning of scratch buffer */
    lc3_scratch_t new_scratch = (lc3_scratch_t) start_address;
    new_scratch->max_scratch_calculation_only = max_scratch_calculation;
    new_scratch->current_alloc_size = 0;
    new_scratch->max_alloc_size = 0;

#ifdef VERBOSE_SCRATCH_ALLOC
    new_scratch->current_alloc_size = offset;
    new_scratch->max_address = (uintptr_t) start + size - 1;
    new_scratch->stack_index = 0;
    new_scratch->max_stack_index = 0;
    new_scratch->alloc_sizes[0] = 0;
#endif

    offset = sizeof( *new_scratch );
    if ( offset & ( alignment - 1 ) )
    {
        offset += alignment - ( offset & ( alignment - 1 ) );
    }

    if ( size < offset )
    {
#ifdef VERBOSE_SCRATCH_ALLOC
        fprintf( stderr, "warning: scratch size too small in lc3_scratch_init" );
        assert( 0 );
#endif
        return NULL;
    }

    if ( max_scratch_calculation
    )
    {
        new_scratch->current_alloc_size += offset;
        new_scratch->max_alloc_size = new_scratch->current_alloc_size;
    }

    new_scratch->current_block_address = start_address + offset;
    new_scratch->current_block_size = 0;
    new_scratch->alignment = alignment;

    size_t address_size = sizeof( new_scratch->current_block_size );
    address_size = ( address_size + alignment - 1 ) & ~( alignment - 1 );

    new_scratch->address_size = address_size;

    return new_scratch;
}

#if defined( VERBOSE_SCRATCH_ALLOC ) || defined( DEBUG_MEM )
void* _lc3_scratch_push( lc3_scratch_t scratch, size_t size, const char* func, const char* file, int line )
#else
void* lc3_scratch_push( lc3_scratch_t scratch, size_t size )
#endif
{
#ifdef VERBOSE_SCRATCH_ALLOC
    size_t size_aligned = ( size + 4 + scratch->alignment - 1 ) & ~( scratch->alignment - 1 );
#else
    size_t size_aligned = ( size + scratch->alignment - 1 ) & ~( scratch->alignment - 1 );
#  ifdef DEBUG_MEM
    UNUSED( func );
    UNUSED( file );
    UNUSED( line );
#  endif
#endif
    uintptr_t new_address = scratch->current_block_address + scratch->current_block_size + scratch->address_size;
    lc3_block_size_t* old_size_ptr = (lc3_block_size_t*) ( scratch->current_block_address + scratch->current_block_size );

    /* push last size on stack */
    *old_size_ptr = scratch->current_block_size;
    scratch->current_block_size = (lc3_block_size_t) size_aligned;
    scratch->current_block_address = new_address;

    if ( scratch->max_scratch_calculation_only
    )
    {
        scratch->current_alloc_size += scratch->current_block_size + scratch->address_size;
        scratch->max_alloc_size = MAX( scratch->current_alloc_size, scratch->max_alloc_size );
    }

#ifdef VERBOSE_SCRATCH_ALLOC
    if ( scratch->current_block_address + scratch->current_block_size - 1 > scratch->max_address )
    {
        fprintf( stderr, "warning: maximal scratch size exceeded in scratch_push; called from %s:%s:%d\n", func, file, line );
        fprintf( stderr, "        tried to allocate %u,   only  %d available \n", (unsigned int) scratch->current_block_size, (int) ( scratch->max_address - scratch->current_block_address ) );
        assert( 0 );
        return NULL;
    }

    char* msg = (char*) new_address + size_aligned - 4;
    strcpy( msg, "LC3" );
    scratch->stack_index++;
    scratch->max_stack_index = MAX( scratch->stack_index, scratch->max_stack_index );
    scratch->alloc_sizes[scratch->stack_index] = (lc3_block_size_t) size_aligned;

#  ifdef VERBOSE_SCRATCH_ALLOC
    fprintf( stderr, "pushing    %-5zu bytes at address %" PRIuPTR "; current_alloc_size = %-5zu; called from %s:%s:%d\n", size_aligned, new_address, scratch->current_alloc_size, func, file, line );
#  endif
#endif /* DEBUG */

    return (void*) new_address;
}

#if defined( VERBOSE_SCRATCH_ALLOC ) || defined( DEBUG_MEM )
void* _lc3_scratch_pop( lc3_scratch_t scratch, void* ptr, const char* func, const char* file, int line )
#else
void* lc3_scratch_pop( lc3_scratch_t scratch, void* ptr )
#endif
{
    lc3_block_size_t* new_size_ptr = (lc3_block_size_t*) ( scratch->current_block_address - scratch->address_size );

    if ( scratch->max_scratch_calculation_only
    )
    {
        scratch->current_alloc_size -= scratch->current_block_size + scratch->address_size;
    }

#ifdef VERBOSE_SCRATCH_ALLOC
    uintptr_t address = (uintptr_t) ptr;

    if ( address != scratch->current_block_address )
    {
        fprintf( stderr, "warning: invalid address %" PRIuPTR " in scratch_pop; expected %" PRIuPTR "; called from %s:%s:%d\n", address, scratch->current_block_address, func, file, line );
        assert( 0 );
        return (void*) address;
    }

    if ( scratch->stack_index == 0 )
    {
        fprintf( stderr, "warning: call to lc3_scratch_pop without previous call to lc3_scratch_push; called from %s:%s:%d\n", func, file, line );
        assert( 0 );
        return (void*) address;
    }

    /* check message integrity */
    const char* msg = (const char*) address + scratch->current_block_size - 4;
    if ( strcmp( "LC3", msg ) )
    {
        fprintf( stderr, "memory access violation detected in region starting with address %" PRIuPTR "; called from %s:%s:%d (msg: %c%c%c%c)\n", address, func, file, line, msg[0], msg[1], msg[2], msg[3] );
        assert( 0 );
    }

#  ifdef VERBOSE_SCRATCH_ALLOC
    fprintf( stderr, "popping    %-5u bytes at address %" PRIuPTR "; current_alloc_size = %-5zu; called from %s:%s:%d\n", scratch->current_block_size,
             address, scratch->current_alloc_size, func, file, line );
#  endif

    scratch->stack_index--;

    if ( *new_size_ptr != scratch->alloc_sizes[scratch->stack_index] )
    {
        fprintf( stderr, "warning: current_block_size changed on stack, check for underflows in this buffer or overflow in previous buffers; called from %s:%s:%d.\n", func, file, line );
        assert( 0 );
        *new_size_ptr = scratch->alloc_sizes[scratch->stack_index];  // should this be done?
    }
#else
#  ifdef DEBUG_MEM
    UNUSED( func );
    UNUSED( file );
    UNUSED( line );
#  endif
    UNUSED( ptr );
#endif

    scratch->current_block_size = *new_size_ptr;
    scratch->current_block_address -= scratch->current_block_size + scratch->address_size;
    return NULL;
}

#ifdef VERBOSE_SCRATCH_ALLOC
void lc3_scratch_print_status( lc3_scratch_t scratch )
{
    fprintf( stderr, "allocated blocks     : %d\n", scratch->stack_index );
    fprintf( stderr, "max size used        : %-5zu\n", scratch->max_alloc_size );
    fprintf( stderr, "max blocks allocated : %-5d\n", scratch->max_stack_index );
}
#endif
