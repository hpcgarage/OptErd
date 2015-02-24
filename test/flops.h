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

#pragma once

#include <config.h>

#if ERD_RECORD_FLOPS
    #include <inttypes.h>
    #include <errno.h>
    #include <unistd.h>
    #include <sys/ioctl.h>
    #include <linux/perf_event.h>
    #include <asm/unistd.h>

    static int perf_event_open(struct perf_event_attr *hw_event, pid_t pid, int cpu, int group_fd, unsigned long flags) {
        return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
    }

    volatile uint64_t flops = 0;
#endif



#if ERD_RECORD_FLOPS
    #if ERD_RECORD_K15_DPADD_FLOPS
        #define PERF_ATTR_CONFIG 0x1003
        #define PERF_COUNTER_NAME "DP ADD GFLOPs"
    #elif ERD_RECORD_K15_DPMUL_FLOPS
        #define PERF_ATTR_CONFIG 0x2003
        #define PERF_COUNTER_NAME "DP MUL GFLOPs"
    #elif ERD_RECORD_K15_DPFMA_FLOPS
        #define PERF_ATTR_CONFIG 0x8003
        #define PERF_COUNTER_NAME "DP FMA GFLOPs"
    #elif ERD_RECORD_K15_DPDIVSQRT_FLOPS
        #define PERF_ATTR_CONFIG 0x4003
        #define PERF_COUNTER_NAME "DP DIV/SQRT GFLOPs"
    #elif ERD_RECORD_K15_DPANY_FLOPS
        #define PERF_ATTR_CONFIG 0xF003
        #define PERF_COUNTER_NAME "DP GFLOPs"
    #elif ERD_RECORD_C2D_PACKEDDP_UOPS
        #define PERF_ATTR_CONFIG 0x04CA
        #define PERF_COUNTER_NAME "Packed DP GuOPs"
    #elif ERD_RECORD_C2D_SCALARDP_UOPS
        #define PERF_ATTR_CONFIG 0x08CA
        #define PERF_COUNTER_NAME "Scalar DP GuOPs"
    #elif ERD_RECORD_WSM_ISSUED_PACKEDFP_UOPS
        #define PERF_ATTR_CONFIG 0x1010
        #define PERF_COUNTER_NAME "Issued Packed FP GuOPs"
    #elif ERD_RECORD_WSM_ISSUED_SCALARFP_UOPS
        #define PERF_ATTR_CONFIG 0x2010
        #define PERF_COUNTER_NAME "Issued Scalar FP GuOPs"
    #elif ERD_RECORD_WSM_PACKEDDP_UOPS
        #define PERF_ATTR_CONFIG 0x04C7
        #define PERF_COUNTER_NAME "Packed DP GuOPs"
    #elif ERD_RECORD_WSM_SCALARDP_UOPS
        #define PERF_ATTR_CONFIG 0x08C7
        #define PERF_COUNTER_NAME "Scalar DP GuOPs"
    #else
        #error "Implement me"
    #endif

    #define BEGIN_RECORD_FLOPS \
        struct perf_event_attr flops_attr; \
        memset(&flops_attr, 0, sizeof(struct perf_event_attr)); \
        flops_attr.size = sizeof(struct perf_event_attr); \
        flops_attr.type = PERF_TYPE_RAW; \
        flops_attr.config = PERF_ATTR_CONFIG; \
        flops_attr.disabled = 1; \
        flops_attr.exclude_kernel = 1; \
        flops_attr.exclude_hv = 1; \
        const int flops_fd = perf_event_open(&flops_attr, 0, -1, -1, 0); \
        if (flops_fd == -1) { \
            fprintf(stderr, "Failed to create perf counter (errno = %d)\n", errno); \
        } \
        assert(ioctl(flops_fd, PERF_EVENT_IOC_RESET, 0) == 0); \
        assert(ioctl(flops_fd, PERF_EVENT_IOC_ENABLE, 0) == 0);
    #define END_RECORD_FLOPS \
        { \
            assert(ioctl(flops_fd, PERF_EVENT_IOC_DISABLE, 0) == 0); \
            uint64_t counter; \
            assert(read(flops_fd, &counter, sizeof(uint64_t)) == sizeof(uint64_t)); \
            __sync_fetch_and_add(&flops, counter); \
            assert(close(flops_fd) == 0); \
        }
    #define REPORT_FLOPS \
        printf("%s: %.3lf\n", PERF_COUNTER_NAME, ((double)flops) * 1.0e-9);
#else
    #define BEGIN_RECORD_FLOPS
    #define END_RECORD_FLOPS
    #define REPORT_FLOPS
#endif
