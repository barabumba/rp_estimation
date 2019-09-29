/*
 * misc.c
 *
 *  Created on: Sep 29, 2019
 *      Author: mturbin
 */
#include "defs.h"


uint64_t GetRandomSeed_64b(void)
{
	uint64_t data;
	FILE *fp;
	fp = fopen("/dev/urandom", "r");
	fread(&data, sizeof(data), 1, fp);
	fclose(fp);
	return data;
}
