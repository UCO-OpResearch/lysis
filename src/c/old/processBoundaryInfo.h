#pragma once
#include "all.h"

typedef struct {
	int processNumber;
	bool bottomBoundary;
	bool leftBoundary;
	bool topBoundary;
	bool rightBoundary;
} ProcessBoundaryInfo;

void setGridBoundaries(ProcessBoundaryInfo* pbInfo, RankGrid* procGrid);
void outputGridBoundaries(ProcessBoundaryInfo* pbInfo);
