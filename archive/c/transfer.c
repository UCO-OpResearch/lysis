
#include "all.h"

void zeroArray(int* array, int count) {
	for (int iii = 0; iii < count; iii++) {
		array[iii] = 0;
	}
}

void bufferIntoFibers_Up(NodeGrid* grid, Parameters* parameters, int time, int* buffer) {
	int bufferIndex = 0;
	int nodeIndex = grid->indexOfLastGhostNode + 1;
	while (bufferIndex < grid->transferUnboundToUpArrayCount) {
		//grid->nodeArray[nodeIndex].upFiber.unboundMolecules += buffer[bufferIndex];
		addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].upFiber, parameters, time, buffer[bufferIndex]);
		bufferIndex++;
		//grid->nodeArray[nodeIndex].rightFiber.unboundMolecules += buffer[bufferIndex];
		addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].rightFiber, parameters, time, buffer[bufferIndex]);
		bufferIndex++;
		nodeIndex++;
	}
}

void bufferIntoFibers_Down(NodeGrid* grid, Parameters* parameters, int time, int* buffer) {
	int bufferIndex = 0;
	//int nodeIndex = grid->nodeWidth * (grid->nodeLengthWithGhostRows - 1);
	int nodeIndex = grid->totalNodesWithGhostRows - grid->nodeWidth;
	while (bufferIndex < grid->transferUnboundToDownArrayCount) {
		//grid->nodeArray[nodeIndex].forwardFiber.unboundMolecules += buffer[bufferIndex];
		addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].forwardFiber, parameters, time, buffer[bufferIndex]);
		bufferIndex++;
		nodeIndex++;
	}
}

void bufferIntoFibers_Right(NodeGrid* grid, Parameters* parameters, int time, int* buffer) {
	int bufferIndex = 0;
	int nodeIndex = 0;
	while (bufferIndex < grid->transferUnboundToRightArrayCount) {
		//grid->nodeArray[nodeIndex].forwardFiber.unboundMolecules += buffer[bufferIndex];
		addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].forwardFiber, parameters, time, buffer[bufferIndex]);
		bufferIndex++;
		//grid->nodeArray[nodeIndex].upFiber.unboundMolecules += buffer[bufferIndex];
		addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].upFiber, parameters, time, buffer[bufferIndex]);
		bufferIndex++;
		nodeIndex += grid->nodeWidth;
	}
}

void bufferIntoFibers_Left(NodeGrid* grid, Parameters* parameters, int time, int* buffer) {
	int bufferIndex = 0;
	int nodeIndex = grid->nodeWidth - 1;
	while (bufferIndex < grid->transferUnboundToLeftArrayCount) {
		//grid->nodeArray[nodeIndex].rightFiber.unboundMolecules += buffer[bufferIndex];
		addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].rightFiber, parameters, time, buffer[bufferIndex]);
		bufferIndex++;
		nodeIndex += grid->nodeWidth;
	}
}

void transferToOtherRanks(NodeGrid* grid, RankGrid* rankGrid, Parameters* parameters, int time) {

	/*****************************/
	// up transfers

	if (rankGrid->length > 1) {

		if (!grid->topBoundary && EVEN == grid->evenOrOddRow) {
			MPI_Send(grid->transferUnboundToUpArray, grid->transferUnboundToUpArrayCount, MPI_INT, grid->MPIRank + rankGrid->width, 0, MPI_COMM_WORLD);
		}

		if (!grid->bottomBoundary && ODD == grid->evenOrOddRow) {
			int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToUpArrayCount);
			MPI_Recv(buffer, grid->transferUnboundToUpArrayCount, MPI_INT, grid->MPIRank - rankGrid->width, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			bufferIntoFibers_Up(grid, parameters, time, buffer);
			free(buffer);
		}

		if (!grid->topBoundary && ODD == grid->evenOrOddRow) {
			MPI_Send(grid->transferUnboundToUpArray, grid->transferUnboundToUpArrayCount, MPI_INT, grid->MPIRank + rankGrid->width, 0, MPI_COMM_WORLD);
		}

		if (!grid->bottomBoundary && EVEN == grid->evenOrOddRow) {
			int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToUpArrayCount);
			MPI_Recv(buffer, grid->transferUnboundToUpArrayCount, MPI_INT, grid->MPIRank - rankGrid->width, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			bufferIntoFibers_Up(grid, parameters, time, buffer);
			free(buffer);
		}

		/*****************************/
		//down transfers

		if (!grid->bottomBoundary && EVEN == grid->evenOrOddRow) {
			MPI_Send(grid->transferUnboundToDownArray, grid->transferUnboundToDownArrayCount, MPI_INT, grid->MPIRank - rankGrid->width, 0, MPI_COMM_WORLD);
		}

		if (!grid->topBoundary && ODD == grid->evenOrOddRow) {
			int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToDownArrayCount);
			MPI_Recv(buffer, grid->transferUnboundToDownArrayCount, MPI_INT, grid->MPIRank + rankGrid->width, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			bufferIntoFibers_Down(grid, parameters, time, buffer);
			free(buffer);
		}

		if (!grid->bottomBoundary && ODD == grid->evenOrOddRow) {
			MPI_Send(grid->transferUnboundToDownArray, grid->transferUnboundToDownArrayCount, MPI_INT, grid->MPIRank - rankGrid->width, 0, MPI_COMM_WORLD);
		}

		if (!grid->topBoundary && EVEN == grid->evenOrOddRow) {
			int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToDownArrayCount);
			MPI_Recv(buffer, grid->transferUnboundToDownArrayCount, MPI_INT, grid->MPIRank + rankGrid->width, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			bufferIntoFibers_Down(grid, parameters, time, buffer);
			free(buffer);
		}

	}

	/*****************************/
	//right transfers

	if (!grid->rightBoundary && 0 == grid->MPIRank % 2) {
		MPI_Send(grid->transferUnboundToRightArray, grid->transferUnboundToRightArrayCount, MPI_INT, grid->MPIRank + 1, 0, MPI_COMM_WORLD);
	}

	if (!grid->leftBoundary && 0 != grid->MPIRank % 2) {
		int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToRightArrayCount);
		MPI_Recv(buffer, grid->transferUnboundToRightArrayCount, MPI_INT, grid->MPIRank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		bufferIntoFibers_Right(grid, parameters, time, buffer);
		free(buffer);
	}

	if (!grid->rightBoundary && 0 != grid->MPIRank % 2) {
		MPI_Send(grid->transferUnboundToRightArray, grid->transferUnboundToRightArrayCount, MPI_INT, grid->MPIRank + 1, 0, MPI_COMM_WORLD);
	}

	if (!grid->leftBoundary && 0 == grid->MPIRank % 2) {
		int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToRightArrayCount);
		MPI_Recv(buffer, grid->transferUnboundToRightArrayCount, MPI_INT, grid->MPIRank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		bufferIntoFibers_Right(grid, parameters, time, buffer);
		free(buffer);
	}

	/*****************************/
	//left transfers

	if (!grid->leftBoundary && 0 == grid->MPIRank % 2) {
		MPI_Send(grid->transferUnboundToLeftArray, grid->transferUnboundToLeftArrayCount, MPI_INT, grid->MPIRank - 1, 0, MPI_COMM_WORLD);
	}

	if (!grid->rightBoundary && 0 != grid->MPIRank % 2) {
		int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToLeftArrayCount);
		MPI_Recv(buffer, grid->transferUnboundToLeftArrayCount, MPI_INT, grid->MPIRank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		bufferIntoFibers_Left(grid, parameters, time, buffer);
		free(buffer);
	}

	if (!grid->leftBoundary && 0 != grid->MPIRank % 2) {
		MPI_Send(grid->transferUnboundToLeftArray, grid->transferUnboundToLeftArrayCount, MPI_INT, grid->MPIRank - 1, 0, MPI_COMM_WORLD);
	}

	if (!grid->rightBoundary && 0 == grid->MPIRank % 2) {
		int* buffer = (int*)malloc(sizeof(int) * grid->transferUnboundToLeftArrayCount);
		MPI_Recv(buffer, grid->transferUnboundToLeftArrayCount, MPI_INT, grid->MPIRank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		bufferIntoFibers_Left(grid, parameters, time, buffer);
		free(buffer);
	}

	/*****************************/
	//diagonal transfers to upper left

	if (rankGrid->length > 1) {

		if (!grid->leftBoundary && !grid->topBoundary && EVEN == grid->evenOrOddRow) {
			//printf("Rank %i sending to %i\n", grid->MPIRank, grid->MPIRank + rankGrid->width - 1);
			MPI_Send(&grid->transferToUpperLeft, 1, MPI_INT, grid->MPIRank + rankGrid->width - 1, 0, MPI_COMM_WORLD);
		}

		if (!grid->rightBoundary && !grid->bottomBoundary && ODD == grid->evenOrOddRow) {
			int buffer;
			MPI_Recv(&buffer, 1, MPI_INT, grid->MPIRank - rankGrid->width + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//grid->nodeArray[grid->nodeWidth - 1].rightFiber.unboundMolecules += buffer;
			addMultipleUnboundToFiber(&grid->nodeArray[grid->nodeWidth - 1].rightFiber, parameters, time, buffer);
		}

		if (!grid->leftBoundary && !grid->topBoundary && ODD == grid->evenOrOddRow) {
			//printf("Rank %i sending to %i\n", grid->MPIRank, grid->MPIRank + rankGrid->width - 1);
			MPI_Send(&grid->transferToUpperLeft, 1, MPI_INT, grid->MPIRank + rankGrid->width - 1, 0, MPI_COMM_WORLD);
		}

		if (!grid->rightBoundary && !grid->bottomBoundary && EVEN == grid->evenOrOddRow) {
			int buffer;
			MPI_Recv(&buffer, 1, MPI_INT, grid->MPIRank - rankGrid->width + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//grid->nodeArray[grid->nodeWidth - 1].rightFiber.unboundMolecules += buffer;
			addMultipleUnboundToFiber(&grid->nodeArray[grid->nodeWidth - 1].rightFiber, parameters, time, buffer);
		}

		/*****************************/
		//diagonal transfers to lower right

		if (!grid->rightBoundary && !grid->bottomBoundary && EVEN == grid->evenOrOddRow) {
			//printf("Rank %i sending to %i\n", grid->MPIRank, grid->MPIRank - rankGrid->width + 1);
			MPI_Send(&grid->transferToLowerRight, 1, MPI_INT, grid->MPIRank - rankGrid->width + 1, 0, MPI_COMM_WORLD);
		}

		if (!grid->leftBoundary && !grid->topBoundary && ODD == grid->evenOrOddRow) {
			int buffer;
			MPI_Recv(&buffer, 1, MPI_INT, grid->MPIRank + rankGrid->width - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			int nodeIndex = grid->nodeWidth * (grid->nodeLengthWithGhostRows - 1);
			//grid->nodeArray[nodeIndex].forwardFiber.unboundMolecules += buffer;
			addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].forwardFiber, parameters, time, buffer);
		}

		if (!grid->rightBoundary && !grid->bottomBoundary && ODD == grid->evenOrOddRow) {
			//printf("Rank %i sending to %i\n", grid->MPIRank, grid->MPIRank - rankGrid->width + 1);
			MPI_Send(&grid->transferToLowerRight, 1, MPI_INT, grid->MPIRank - rankGrid->width + 1, 0, MPI_COMM_WORLD);
		}

		if (!grid->leftBoundary && !grid->topBoundary && EVEN == grid->evenOrOddRow) {
			int buffer;
			MPI_Recv(&buffer, 1, MPI_INT, grid->MPIRank + rankGrid->width - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			int nodeIndex = grid->nodeWidth * (grid->nodeLengthWithGhostRows - 1);
			//grid->nodeArray[nodeIndex].forwardFiber.unboundMolecules += buffer;
			addMultipleUnboundToFiber(&grid->nodeArray[nodeIndex].forwardFiber, parameters, time, buffer);
		}
	}

	/*****************************/
	// reset the transfer arrays

	if (NULL != grid->transferUnboundToUpArray) {
		zeroArray(grid->transferUnboundToUpArray, grid->transferUnboundToUpArrayCount);
	}
	if (NULL != grid->transferUnboundToRightArray) {
		zeroArray(grid->transferUnboundToRightArray, grid->transferUnboundToRightArrayCount);
	}
	if (NULL != grid->transferUnboundToDownArray) {
		zeroArray(grid->transferUnboundToDownArray, grid->transferUnboundToDownArrayCount);
	}
	if (NULL != grid->transferUnboundToLeftArray) {
		zeroArray(grid->transferUnboundToLeftArray, grid->transferUnboundToLeftArrayCount);
	}

	grid->transferToLowerRight = 0;
	grid->transferToUpperLeft = 0;
}

/*****************************************************************/
/*****************************************************************/

void transferInputData(InputData* inputData, Parameters* parameters) {
	MPI_Bcast(inputData->unbindingTimeDistribution, parameters->unbindingTimeDist, MPI_DOUBLE, HEAD_RANK, MPI_COMM_WORLD);
	for (int index = 0; index < parameters->maxLysesPerBlock; index++) {
		MPI_Bcast(inputData->lysisTime[index], parameters->lysisBlocks, MPI_DOUBLE, HEAD_RANK, MPI_COMM_WORLD);
	}
	MPI_Bcast(inputData->lysesPerBlock, parameters->lysisBlocks, MPI_INT, HEAD_RANK, MPI_COMM_WORLD);
 }

/*****************************************************************/
/*****************************************************************/

/*void test_MarkTransferArrays(NodeGrid* grid) {

	int mark = grid->MPIRank * 1000;

	if (NULL != grid->transferUnboundToUpArray) {
		for (int iii = 0; iii < grid->transferUnboundToUpArrayCount; iii++) {
			grid->transferUnboundToUpArray[iii] = mark;
		}
	}
	if (NULL != grid->transferUnboundToRightArray) {
		for (int iii = 0; iii < grid->transferUnboundToRightArrayCount; iii++) {
			grid->transferUnboundToRightArray[iii] = mark;
		}
	}
	if (NULL != grid->transferUnboundToDownArray) {
		for (int iii = 0; iii < grid->transferUnboundToDownArrayCount; iii++) {
			grid->transferUnboundToDownArray[iii] = mark;
		}
	}
	if (NULL != grid->transferUnboundToLeftArray) {
		for (int iii = 0; iii < grid->transferUnboundToLeftArrayCount; iii++) {
			grid->transferUnboundToLeftArray[iii] = mark;
		}
	}
	if (!grid->rightBoundary && !grid->bottomBoundary) {
		grid->transferToLowerRight = mark;
	}

	if (!grid->leftBoundary && !grid->topBoundary) {
		grid->transferToUpperLeft = mark;
	}
}*/