
#include "all.h"

void retrieveFromFibers(NodeGrid* grid, RankGrid* rankGrid, int* localBuffer, char* fileName);

void printUnboundMoleculesToFile(NodeGrid* grid, RankGrid* rankGrid, char* fileName) {

	int* localBuffer = NULL;
	localBuffer = (int*)malloc(sizeof(int) * grid->totalFibersWithGhostRows);

	int fiber = 0;
	for (int iii = 0; iii < grid->totalNodesWithGhostRows; iii++) {
		localBuffer[fiber] = grid->nodeArray[iii].forwardFiber.unboundMolecules;
		fiber++;
		localBuffer[fiber] = grid->nodeArray[iii].upFiber.unboundMolecules;
		fiber++;
		localBuffer[fiber] = grid->nodeArray[iii].rightFiber.unboundMolecules;
		fiber++;
	}

	retrieveFromFibers(grid, rankGrid, localBuffer, fileName);
	free(localBuffer);
}

void printBoundMoleculesToFile(NodeGrid* grid, RankGrid* rankGrid, char* fileName) {

	int* localBuffer = NULL;
	localBuffer = (int*)malloc(sizeof(int) * grid->totalFibersWithGhostRows);

	int fiber = 0;
	for (int iii = 0; iii < grid->totalNodesWithGhostRows; iii++) {
		localBuffer[fiber] = grid->nodeArray[iii].forwardFiber.boundMolecules;
		fiber++;
		localBuffer[fiber] = grid->nodeArray[iii].upFiber.boundMolecules;
		fiber++;
		localBuffer[fiber] = grid->nodeArray[iii].rightFiber.boundMolecules;
		fiber++;
	}

	retrieveFromFibers(grid, rankGrid, localBuffer, fileName);
	free(localBuffer);
}

void printIsDegradedToFile(NodeGrid* grid, RankGrid* rankGrid, char* fileName) {

	int* localBuffer = NULL;
	localBuffer = (int*)malloc(sizeof(int) * grid->totalFibersWithGhostRows);

	int fiber = 0;
	for (int iii = 0; iii < grid->totalNodesWithGhostRows; iii++) {
		localBuffer[fiber] = grid->nodeArray[iii].forwardFiber.fiberIsDegraded;
		fiber++;
		localBuffer[fiber] = grid->nodeArray[iii].upFiber.fiberIsDegraded;
		fiber++;
		localBuffer[fiber] = grid->nodeArray[iii].rightFiber.fiberIsDegraded;
		fiber++;
	}

	retrieveFromFibers(grid, rankGrid, localBuffer, fileName);
	free(localBuffer);
}



void retrieveFromFibers(NodeGrid* grid, RankGrid* rankGrid, int* localBuffer, char* fileName) {

	int rank;
	int size;
	//MPI_Status status;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/***********************************************************************/
	//send count

	int totalFibers = 0;
	MPI_Reduce(&grid->totalFibersWithGhostRows, &totalFibers, 1, MPI_INT, MPI_SUM, HEAD_RANK, MPI_COMM_WORLD);
	
	int* receiveBuffer = NULL;
	int* receiveTotalFibersPerRank = NULL;
	int* receiveDisplacement = NULL;
	int* receiveRowFiberWidthPerRank = NULL;
	int* receiveRowsPerRank = NULL;
	if (HEAD_RANK == rank) {
		receiveBuffer = (int*)malloc(sizeof(int) * totalFibers);
		receiveTotalFibersPerRank = (int*)malloc(sizeof(int) * size);
		receiveDisplacement = (int*)malloc(sizeof(int) * size);
		receiveRowFiberWidthPerRank = (int*)malloc(sizeof(int) * size);
		receiveRowsPerRank = (int*)malloc(sizeof(int) * size);
	}

	MPI_Gather(&grid->totalFibersWithGhostRows, 1, MPI_INT, receiveTotalFibersPerRank, 1, MPI_INT, HEAD_RANK, MPI_COMM_WORLD);
	MPI_Gather(&grid->fiberWidth, 1, MPI_INT, receiveRowFiberWidthPerRank, 1, MPI_INT, HEAD_RANK, MPI_COMM_WORLD);
	MPI_Gather(&grid->nodeLengthWithGhostRows, 1, MPI_INT, receiveRowsPerRank, 1, MPI_INT, HEAD_RANK, MPI_COMM_WORLD);

	//calculate array position for each rank
	if (HEAD_RANK == rank) {
		receiveDisplacement[0] = 0;
		for (int iii = 1; iii < size; iii++) {
			receiveDisplacement[iii] = receiveTotalFibersPerRank[iii - 1] + receiveDisplacement[iii - 1];
		}
	}

	MPI_Gatherv(localBuffer, grid->totalFibersWithGhostRows, MPI_INT, receiveBuffer, receiveTotalFibersPerRank, receiveDisplacement, MPI_INT, HEAD_RANK, MPI_COMM_WORLD);


	/***********************************************************************/
	//HEAD_RANK print to file

	if (HEAD_RANK == rank) {

		FILE* file = fopen(fileName, "a"); // "w" for new file, "a" for append to existing

		int displacementToRow = 0;
		int startReceiveBufferPosition;
		int endReceiveBufferPosition;

		for (int firstRankInRow = 0; firstRankInRow < rankGrid->total; firstRankInRow += rankGrid->width) {
			for (int row = 0; row < receiveRowsPerRank[firstRankInRow]; row++) {
				for (int recvRank = firstRankInRow; recvRank < firstRankInRow + rankGrid->width; recvRank++) {
					displacementToRow = row * receiveRowFiberWidthPerRank[recvRank];
					startReceiveBufferPosition = displacementToRow + receiveDisplacement[recvRank];
					endReceiveBufferPosition = startReceiveBufferPosition + receiveRowFiberWidthPerRank[recvRank];
					for (int fiber = startReceiveBufferPosition; fiber < endReceiveBufferPosition; fiber++) {
						fprintf(file, "%i,", receiveBuffer[fiber]);
					}
				}
				fprintf(file, "\r\n");
			}
		}
	}

	free(receiveBuffer);
	free(receiveTotalFibersPerRank);
	free(receiveDisplacement);
	free(receiveRowFiberWidthPerRank);
	free(receiveRowsPerRank);
}

//for debugging
/*void test_MarkForInspection(NodeGrid* grid) {
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int mark = rank; //* 1000;

	for (int iii = 0; iii < grid->totalNodesWithGhostRows; iii++) {
		grid->nodeArray[iii].forwardFiber.unboundMolecules = mark;
		//mark++;
		grid->nodeArray[iii].upFiber.unboundMolecules = mark;
		//mark++;
		grid->nodeArray[iii].rightFiber.unboundMolecules = mark;
		//mark++;
	}
}*/