/*
 *
 * $Id: rand.h,v 1.1 2005/01/12 10:07:51 fringer Exp fringer $
 * $Log: rand.h,v $
 * Revision 1.1  2005/01/12 10:07:51  fringer
 * Initial revision
 *
 *
 */
#include<stdlib.h>

#define frand(xmin,xmax) ((double)xmin+(double)(xmax-xmin)*rand()/(double)RAND_MAX)
