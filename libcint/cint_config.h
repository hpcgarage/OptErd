/*
 * Copyright (c) 2013-2015 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * The GNU Lesser General Public License is included in this distribution
 * in the file COPYING.
 */

#ifndef __CINT_CONFIG_H__
#define __CINT_CONFIG_H__


#include <stdio.h>


#ifndef __APPLE__
#define HAS_MALLOC_H
#endif

#define _DEBUG_LEVEL_    1   //  0 to 10,
                             //  0 is no debug print info at all,
                             //  10 is full info

#define ALIGNED_MALLOC(size)  _mm_malloc(size, __ALIGNLEN__)
#define ALIGNED_FREE(addr)    _mm_free(addr)


#if ( _DEBUG_LEVEL_ == -1 )
#define CINT_PRINTF( level, fmt, args... )        {}
#else
#define CINT_PRINTF( level, fmt, args... )                                          \
        do                                                                          \
        {                                                                           \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ )                               \
            {                                                                       \
                sprintf( basis->str_buf, "%s() line %d ", __FUNCTION__, __LINE__ ); \
                sprintf( basis->str_buf + strlen(basis->str_buf), fmt, ##args );    \
                fprintf( stdout, "%s", basis->str_buf );                            \
                fflush( stdout );                                                   \
            }                                                                       \
        } while ( 0 )
#endif


#if ( _DEBUG_LEVEL_ > 1 )
#define CINT_INFO( fmt, args... )                                            \
        do                                                                   \
        {                                                                    \
            sprintf( basis->str_buf, "**** CInt: ");                         \
            sprintf( basis->str_buf + strlen("**** CInt: "), fmt, ##args );  \
            fprintf( stdout, "%s", basis->str_buf );                         \
            fflush( stdout );                                                \
        } while ( 0 )
#else
#define CINT_INFO( fmt, args... )        {}
#endif


#endif /* __CINT_CONFIG_H__ */
