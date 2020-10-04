/*
 * defs.h
 *
 *  Created on: May 1, 2019
 *      Author: mturbin
 */

#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 0x4000
#define M 0x800
#define Q 10.

#define M_SRCH_RANGE_LEN 0x80
#define Q_SRCH_RANGE_LEN 0x40
#define M_SRCH_MIN M/2
#define M_SRCH_MAX 3*M/2
#define Q_SRCH_MIN Q/2
#define Q_SRCH_MAX 3*Q/2

#define NUMBER_OF_TESTS 1000000 

#endif /* SRC_DEFS_H_ */
