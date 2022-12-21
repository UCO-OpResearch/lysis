//#include "processBoundaryInfo.h"
//#include "parameters.h"
#include "all.h"

void setGridBoundaries(ProcessBoundaryInfo* pbInfo, RankGrid* procGrid) {

	//check for left or right boundaries
	if (0 == pbInfo->processNumber % procGrid->width) {
		pbInfo->leftBoundary = true;
		pbInfo->rightBoundary = false;
	}
	else if ((procGrid->width - 1) == pbInfo->processNumber % procGrid->width) {
		pbInfo->leftBoundary = false;
		pbInfo->rightBoundary = true;
	}
	else {
		pbInfo->leftBoundary = false;
		pbInfo->rightBoundary = false;
	}

	//check for top or bottom boundaries
	if (pbInfo->processNumber < procGrid->width) {
		pbInfo->bottomBoundary = true;
		pbInfo->topBoundary = false;
	}
	else if(pbInfo->processNumber > (procGrid->total - procGrid->width - 1)){
		pbInfo->bottomBoundary = false;
		pbInfo->topBoundary = true;
	}
	else {
		pbInfo->bottomBoundary = false;
		pbInfo->topBoundary = false;
	}

}

//output to command line for testing
void outputGridBoundaries(ProcessBoundaryInfo* pbInfo) {
	printf("Process number (pbInfo->processinfo): %i\n", pbInfo->processNumber);
	printf("\tbottom boundary: %i\n", pbInfo->bottomBoundary);
	printf("\ttop boundary: %i\n", pbInfo->topBoundary);
	printf("\tleft boundary: %i\n", pbInfo->leftBoundary);
	printf("\tright boundary: %i\n\n", pbInfo->rightBoundary);
}