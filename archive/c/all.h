#pragma once

#define _GNU_SOURCE

#include <stdint.h>	//for UINT_LEAST32_T
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef NEWTIMELIST
#include "dynamicArray.h"
#else
#include "timeList.h"
#endif

#include "initializeData.h"
#include "mpi.h"
#include "nodeGrid.h"
#include "parameters.h"
#include "retrieve.h"
#include "simulation.h"
#include "test.h"
#include "transfer.h"
#include "kiss.h"

#include <execinfo.h>
#include <signal.h>
#include <unistd.h>