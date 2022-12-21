#pragma once

#include "parameters.h"
#include "timeList.h"

#define EVEN 0
#define ODD 1

#define FORWARD 0
#define UP 1
#define RIGHT 2

typedef struct {
	bool fiberIsDegraded;
	int degradeTime;

	int unboundMolecules;
#ifdef NEWTIMELIST
	DynArray* unBoundTimeList;
#else
	TimeList* unBoundTimeList;
#endif

	int boundMolecules;
#ifdef NEWTIMELIST
	DynArray* boundTimeList;
#else
	TimeList* boundTimeList;
#endif
} Fiber;

typedef struct {
	Fiber forwardFiber; //along the length
	Fiber upFiber; //along the height
	Fiber rightFiber; //along the width
} Node;

typedef struct {
	int MPIRank;
	bool bottomBoundary;
	bool leftBoundary;
	bool topBoundary;
	bool rightBoundary;
	int evenOrOddRow;

	int nodeLength; //forward
	int nodeWidth; //left right
	int nodeHeight;
	int totalNodesWithGhostRows;
	int fiberWidth;
	int totalFibersWithGhostRows;
	int ghostRows; //number of ghost rows, only for the bottom boundary processes
	int nodeLengthWithGhostRows;
	int indexOfLastGhostNode;
	Node* nodeArray;


	int* transferUnboundToUpArray;
	int transferUnboundToUpArrayCount;
	int* transferUnboundToRightArray;
	int transferUnboundToRightArrayCount;
	int* transferUnboundToDownArray;
	int transferUnboundToDownArrayCount;
	int* transferUnboundToLeftArray;
	int transferUnboundToLeftArrayCount;
	int transferToUpperLeft;
	int transferToLowerRight;

} NodeGrid;

void setGridBoundaries(NodeGrid* grid, RankGrid* rankGrid);
void initializeNodeGrid(NodeGrid* grid, TotalGridSize* totalGridSize, Parameters* parameters, RankGrid* rankGrid);
void intializeNodeGridSize(TotalGridSize* totalGridSize, RankGrid* rankGrid, NodeGrid* grid);