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

#if ERD_RECORD_RATIO
    #include <inttypes.h>
    #include <errno.h>
    #include <unistd.h>
    #include <sys/ioctl.h>
    #include <asm/unistd.h>

# ifdef __linux__
    #include <linux/perf_event.h>

    static int perf_event_open(struct perf_event_attr *hw_event, pid_t pid, int cpu, int group_fd, unsigned long flags) {
        return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
    }
# endif

    volatile uint64_t numerator = 0;
    volatile uint64_t denominator = 0;
#endif



#if ERD_RECORD_RATIO
    #if ERD_RECORD_CPI_RATIO
        #define NUMERATOR_PERF_ATTR_TYPE PERF_TYPE_HARDWARE
        #define NUMERATOR_PERF_ATTR_CONFIG PERF_COUNT_HW_CPU_CYCLES
        #define NUMERATOR_PERF_COUNTER_NAME "Cycles"
        #define DENOMINATOR_PERF_ATTR_TYPE PERF_TYPE_HARDWARE
        #define DENOMINATOR_PERF_ATTR_CONFIG PERF_COUNT_HW_INSTRUCTIONS
        #define DENOMINATOR_PERF_COUNTER_NAME "Instructions"
        #define RATIO_NAME "Cycles/Instruction"
        #define RATIO_UNITS ""
        #define RATIO_MULTIPLIER 1.0
    #elif ERD_RECORD_BRMISS_RATIO
        #define NUMERATOR_PERF_ATTR_TYPE PERF_TYPE_HARDWARE
        #define NUMERATOR_PERF_ATTR_CONFIG PERF_COUNT_HW_BRANCH_MISSES
        #define NUMERATOR_PERF_COUNTER_NAME "Mispredictions"
        #define DENOMINATOR_PERF_ATTR_TYPE PERF_TYPE_HARDWARE
        #define DENOMINATOR_PERF_ATTR_CONFIG PERF_COUNT_HW_BRANCH_INSTRUCTIONS
        #define DENOMINATOR_PERF_COUNTER_NAME "Branches"
        #define RATIO_NAME "Branch mispredictions"
        #define RATIO_UNITS "%%"
        #define RATIO_MULTIPLIER 100.0
    #elif ERD_RECORD_BRANCH_RATIO
        #define NUMERATOR_PERF_ATTR_TYPE PERF_TYPE_HARDWARE
        #define NUMERATOR_PERF_ATTR_CONFIG PERF_COUNT_HW_BRANCH_INSTRUCTIONS
        #define NUMERATOR_PERF_COUNTER_NAME "Branches"
        #define DENOMINATOR_PERF_ATTR_TYPE PERF_TYPE_HARDWARE
        #define DENOMINATOR_PERF_ATTR_CONFIG PERF_COUNT_HW_INSTRUCTIONS
        #define DENOMINATOR_PERF_COUNTER_NAME "Instructions"
        #define RATIO_NAME "Branch ratio"
        #define RATIO_UNITS "%%"
        #define RATIO_MULTIPLIER 100.0
    #elif ERD_RECORD_GPINSTR_RATIO
        #error "Implement me"
    #else
        #error "Implement me"
    #endif

    #define BEGIN_RECORD_RATIO \
        struct perf_event_attr numerator_attr, denominator_attr; \
        memset(&numerator_attr, 0, sizeof(struct perf_event_attr)); \
        memset(&denominator_attr, 0, sizeof(struct perf_event_attr)); \
        numerator_attr.size = sizeof(struct perf_event_attr); \
        numerator_attr.type = NUMERATOR_PERF_ATTR_TYPE; \
        numerator_attr.config = NUMERATOR_PERF_ATTR_CONFIG; \
        numerator_attr.disabled = 1; \
        numerator_attr.exclude_kernel = 1; \
        numerator_attr.exclude_hv = 1; \
        denominator_attr.size = sizeof(struct perf_event_attr); \
        denominator_attr.type = DENOMINATOR_PERF_ATTR_TYPE; \
        denominator_attr.config = DENOMINATOR_PERF_ATTR_CONFIG; \
        denominator_attr.disabled = 1; \
        denominator_attr.exclude_kernel = 1; \
        denominator_attr.exclude_hv = 1; \
        const int numerator_fd = perf_event_open(&numerator_attr, 0, -1, -1, 0); \
        assert(numerator_fd != -1); \
        const int denominator_fd = perf_event_open(&denominator_attr, 0, -1, numerator_fd, 0); \
        assert(denominator_fd != -1); \
        assert(ioctl(numerator_fd, PERF_EVENT_IOC_RESET, 0) == 0); \
        assert(ioctl(denominator_fd, PERF_EVENT_IOC_RESET, 0) == 0); \
        assert(ioctl(numerator_fd, PERF_EVENT_IOC_ENABLE, 0) == 0); \
        assert(ioctl(denominator_fd, PERF_EVENT_IOC_ENABLE, 0) == 0);
    #define END_RECORD_RATIO \
        { \
            assert(ioctl(numerator_fd, PERF_EVENT_IOC_DISABLE, 0) == 0); \
            assert(ioctl(denominator_fd, PERF_EVENT_IOC_DISABLE, 0) == 0); \
            uint64_t counter; \
            assert(read(numerator_fd, &counter, sizeof(uint64_t)) == sizeof(uint64_t)); \
            __sync_fetch_and_add(&numerator, counter); \
            assert(read(denominator_fd, &counter, sizeof(uint64_t)) == sizeof(uint64_t)); \
            __sync_fetch_and_add(&denominator, counter); \
            assert(close(numerator_fd) == 0); \
            assert(close(denominator_fd) == 0); \
        }
    #define REPORT_RATIO \
        printf("%s: %"PRIu64"\n", NUMERATOR_PERF_COUNTER_NAME, numerator); \
        printf("%s: %"PRIu64"\n", DENOMINATOR_PERF_COUNTER_NAME, denominator); \
        printf("%s: %5.2lf " RATIO_UNITS "\n", RATIO_NAME, ((double)numerator) / ((double)denominator) * RATIO_MULTIPLIER);
#else
    #define BEGIN_RECORD_RATIO
    #define END_RECORD_RATIO
    #define REPORT_RATIO
#endif
