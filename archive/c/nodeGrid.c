#include "all.h"

void initializeFiber(Fiber* x, const int totalTimeSteps) {
	x->fiberIsDegraded = false;
	x->degradeTime = totalTimeSteps + 1;

	x->unboundMolecules = 0;
	x->unBoundTimeList = NULL;

	x->boundMolecules = 0;
	x->boundTimeList = NULL;
}

void initializeAllFibers(NodeGrid* grid, const int totalTimeSteps) {
	int totalSize = grid->totalNodesWithGhostRows;
	for (int iii = 0; iii < totalSize; iii++) {
		initializeFiber(&grid->nodeArray[iii].forwardFiber, totalTimeSteps);
		initializeFiber(&grid->nodeArray[iii].rightFiber, totalTimeSteps);
		initializeFiber(&grid->nodeArray[iii].upFiber, totalTimeSteps);
	}
}

void setGridBoundaries(NodeGrid* grid, RankGrid* rankGrid) {

	//check for left or right boundaries
	if (0 == grid->MPIRank % rankGrid->width) {
		grid->leftBoundary = true;
		grid->rightBoundary = false;
	}
	else if ((rankGrid->width - 1) == grid->MPIRank % rankGrid->width) {
		grid->leftBoundary = false;
		grid->rightBoundary = true;
	}
	else {
		grid->leftBoundary = false;
		grid->rightBoundary = false;
	}

	//check for top or bottom boundaries
	if (grid->MPIRank < rankGrid->width) {
		grid->bottomBoundary = true;
		grid->topBoundary = false;
	}
	else if (grid->MPIRank > (rankGrid->total - rankGrid->width - 1)) {
		grid->bottomBoundary = false;
		grid->topBoundary = true;
	}
	else {
		grid->bottomBoundary = false;
		grid->topBoundary = false;
	}
}

void setGhostNodesToDegraded(NodeGrid* grid) {
	for (int iii = 0; iii < grid->indexOfLastGhostNode + 1; iii++) {
		grid->nodeArray[iii].forwardFiber.fiberIsDegraded = true;
		grid->nodeArray[iii].upFiber.fiberIsDegraded = true;
		grid->nodeArray[iii].rightFiber.fiberIsDegraded = true;
	}
}

//output to command line for testing
void outputGridBoundaries(NodeGrid* grid) {
	printf("Rank number (grid->MPIRank): %i\n", grid->MPIRank);
	printf("\tbottom boundary: %i\n", grid->bottomBoundary);
	printf("\ttop boundary: %i\n", grid->topBoundary);
	printf("\tleft boundary: %i\n", grid->leftBoundary);
	printf("\tright boundary: %i\n\n", grid->rightBoundary);
}

void intializeNodeGridSize(TotalGridSize* totalGridSize, RankGrid* rankGrid, NodeGrid* grid) {

	grid->nodeHeight = 1; // will need to change this if the model is expanded into 3d

	// set node length
	if(0 == totalGridSize->nodeLengthTotal % rankGrid->length){ // divides equally, so each rank gets same length
		grid->nodeLength = totalGridSize->nodeLengthTotal / rankGrid->length;
	}
	else{
		int size = totalGridSize->nodeLengthTotal / rankGrid->length;
		
		// calculate the rank row to see if the rank gets extra length
		if(grid->MPIRank > rankGrid->total - totalGridSize->nodeLengthTotal % rankGrid->length * rankGrid->width - 1){
			grid->nodeLength = size + 1;
		}
		else{
			grid->nodeLength = size;
		}
	}

	// set node width
	if(0 == totalGridSize->nodeWidthTotal % rankGrid->width){ // divides equally, so each rank gets same length
		grid->nodeWidth = totalGridSize->nodeWidthTotal / rankGrid->width;
	}
	else{ // spread out the nodes among most ranks, with 1 rank taking up the slack
		int size = totalGridSize->nodeWidthTotal / rankGrid->width;
		
		if((grid->MPIRank % rankGrid->width) < totalGridSize->nodeWidthTotal % rankGrid->width){
			grid->nodeWidth = size + 1;
		}
		else{
			grid->nodeWidth = size;
		}
	}

	if(grid->bottomBoundary){
		grid->ghostRows = totalGridSize->ghostRows;
	}
}

void initializeNodeGrid(NodeGrid* grid, TotalGridSize* totalGridSize, Parameters* parameters, RankGrid* rankGrid) {

	intializeNodeGridSize(totalGridSize, rankGrid, grid);

	/********************************************/
	//transfer arrays and variables

	//bottom
	//setup bottom boundary first since bottomBoundary determines the total amount of rows. Important for left and right boundary grids
	grid->transferUnboundToDownArrayCount = grid->nodeWidth;
	if (grid->bottomBoundary) {
		grid->transferUnboundToDownArray = NULL;
		grid->indexOfLastGhostNode = grid->ghostRows * grid->nodeWidth - 1; 
	}
	else {
		grid->ghostRows = 0;
		grid->indexOfLastGhostNode = -1;
		grid->transferUnboundToDownArray = (int*)calloc(grid->transferUnboundToDownArrayCount, sizeof(int)); //one fiber per node sticking up from grid below
	}

	grid->nodeLengthWithGhostRows = grid->nodeLength + grid->ghostRows;

	//up
	grid->transferUnboundToUpArrayCount = grid->nodeWidth * 2;
	if (grid->topBoundary) {
		grid->transferUnboundToUpArray = NULL;
	}
	else {
		grid->transferUnboundToUpArray = (int*)calloc(grid->transferUnboundToUpArrayCount, sizeof(int)); //2 fibers above each fiber on this grid
	}

	//right
	grid->transferUnboundToRightArrayCount = grid->nodeLengthWithGhostRows * 2;
	if (grid->rightBoundary) {
		grid->transferUnboundToRightArray = NULL;
	}
	else {
		grid->transferUnboundToRightArray = (int*)calloc(grid->transferUnboundToRightArrayCount, sizeof(int)); //2 fibers to the right of each fiber on this grid
	}


	//left
	grid->transferUnboundToLeftArrayCount = grid->nodeLengthWithGhostRows;
	if (grid->leftBoundary) {
		grid->transferUnboundToLeftArray = NULL;
	}
	else {
		grid->transferUnboundToLeftArray = (int*)calloc(grid->transferUnboundToLeftArrayCount, sizeof(int)); //2 fibers to the right of each fiber on this grid
	}

	grid->transferToLowerRight = 0;
	grid->transferToUpperLeft = 0;

	/********************************************/

	grid->totalNodesWithGhostRows = (grid->nodeHeight) * (grid->nodeLength + grid->ghostRows) * (grid->nodeWidth);
	grid->totalFibersWithGhostRows = grid->totalNodesWithGhostRows * FIBERS_PER_NODE;
	grid->fiberWidth = grid->nodeWidth * FIBERS_PER_NODE;
	grid->nodeArray = (Node*)malloc(sizeof(Node) * grid->totalNodesWithGhostRows);
	initializeAllFibers(grid, parameters->totalTimeSteps);
	setGhostNodesToDegraded(grid);

	//am I an even or odd row in the RankGrid
	//useful when doing transfers between ranks
	int row = 0;
	while (1) {
		if (grid->MPIRank < rankGrid->width * (row + 1)) {
			if (0 == row % 2) {
				grid->evenOrOddRow = EVEN;
			}
			else {
				grid->evenOrOddRow = ODD;
			}
			break;
		}
		row++;
	}

}