#ifndef GLOBALS_H
#define GLOBALS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

using namespace std;

extern unsigned short **I_i, **I_j;
extern float mu, delta, **H;

extern struct timeval	StartTime, EndTime;

struct array_max_t{
	int ind;
	float max;
};

#endif
