/*
 * File: timer.h
 * Description: Header file for functions defined in timer.c
 * ---------------------------------------------------------
 * In order to implement the desired timer, it must be defined.
 * The possible timers are as follows:
 *
 * #define GETTIMEOFDAY  // gettimeofday(2)
 * #define CLOCK         // clock(3)
 * #define TIMES         // times(2)
 * #define GETRUSAGE     // getrusage(2)
 *
 * Note that the only useful timer for this benchmark is
 * gettimeofday(2), which has microsecond resolution.
 *
 * CME212 Assignment 5
 * Amit Kushwaha
 * Stanford University
 *
 */
#ifndef _timer_h
#define _timer_h

#include<stdio.h>
#include<stdlib.h>

#define GETTIMEOFDAY

#ifdef CLOCK
#include<time.h>
#endif

#ifdef TIMES

#include<unistd.h>
#include<sys/times.h>
struct tms tms_time;

#endif

#ifdef GETRUSAGE

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define evaltime(rstruct) (double)rstruct.ru_utime.tv_sec \
                      + (double)rstruct.ru_utime.tv_usec*1e-6 

#endif

#ifdef GETTIMEOFDAY

#include<sys/time.h>
#define evaltime(timeval_time) (double)timeval_time.tv_sec \
                             + (double)timeval_time.tv_usec*1e-6 

#endif

// The timer type should be a double
typedef double timeType;

/*
 * Function: Timer
 * Usage: printf("Time = %f\n",Timer()-t0);
 * ----------------------------------------
 * Returns the time in seconds and uses the timer
 * defined in timer.h.
 *
 */
extern timeType Timer(void);

#endif
