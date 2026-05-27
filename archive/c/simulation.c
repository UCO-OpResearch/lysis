#include "all.h"

#define NEWCODE 1

void stop(int processNum) {
	int i = 0;
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	printf("Hostname: %s Process: %i PID: %i waiting for gdb attach\n", hostname, processNum, getpid());
	fflush(stdout);
	while (0 == i) {
		sleep(5);
	}
}

/*
* Returns a random integer between 0 and N-1 inclusive.
*/
int randomInt(int N) {
	double t = urcw1_();		// Call the method from 'kiss32' and store the resulting double
	return (int)(t*N);	// Cast the double to an unsigned int and return it
}

//************************************************//
//************************************************//

void distributeMolecules(NodeGrid* grid, RankGrid* rankGrid, int totalMolecules) {

	int myFiberCount;

	if (grid->bottomBoundary) {
		myFiberCount = grid->nodeWidth * 2;
	}
	else {
		myFiberCount = 0;
	}

	int* recvFibers;
	if(HEAD_RANK == grid->MPIRank){
		recvFibers = (int*)calloc(rankGrid->total, sizeof(int));
	}

	MPI_Gather(&myFiberCount, 1, MPI_INT, recvFibers, 1, MPI_INT, HEAD_RANK, MPI_COMM_WORLD);

	/***********************************************/
	//head randomly assigns molecules and sends to bottomBoundary nodes

	if(HEAD_RANK == grid->MPIRank) {
		int totalFibers = 0;
		for (int rank = 0; rank < rankGrid->width; rank++) { //ignore the other nodes counts, only need bottom boundary nodes
			totalFibers += recvFibers[rank];
		}

		int* start = (int*)calloc(totalFibers, sizeof(int));
		int placement;
		for (int iii = 0; iii < totalMolecules; iii++) {
			placement = randomInt(totalFibers - 1); //don't include the last right fiber, it is not part of the simulation
			start[placement]++;
		}

		int index = 0;
		for (int rank = 0; rank < rankGrid->width; rank++) {
			MPI_Send(&start[index], recvFibers[rank], MPI_INT, rank, 0, MPI_COMM_WORLD);
			index += recvFibers[rank];
		}

		free(start);
		free(recvFibers);
	}

	/***********************************************/

	if (grid->bottomBoundary) {
		int* buffer = (int*)malloc(myFiberCount * sizeof(int));
		MPI_Recv(buffer, myFiberCount, MPI_INT, HEAD_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		int nodeIndex = 0;
		int fiberIndex = 0;
		while (fiberIndex < myFiberCount) {
			grid->nodeArray[nodeIndex].upFiber.unboundMolecules = buffer[fiberIndex];
			fiberIndex++;
			grid->nodeArray[nodeIndex].rightFiber.unboundMolecules = buffer[fiberIndex];
			fiberIndex++;
			nodeIndex++;
		}

		free(buffer);
	}
}

//************************************************//
//************************************************//
//degradation

void checkForDegrade_Fiber(Fiber* fiber, int time) {
	if (!fiber->fiberIsDegraded && fiber->degradeTime <= time) {
		fiber->fiberIsDegraded = true;
		fiber->unboundMolecules = fiber->unboundMolecules + fiber->boundMolecules;
		fiber->boundMolecules = 0;
		if (NULL != fiber->unBoundTimeList) {
#ifdef NEWTIMELIST
			delEntireDA(&fiber->boundTimeList);			
#else
			deleteTimeList(&fiber->boundTimeList);
#endif
		}
	}
}

void checkForDegrade_Node(Node* node, int time) { //checks to see if fibers need to be degraded and degrades if needed
	checkForDegrade_Fiber(&node->forwardFiber, time);
	checkForDegrade_Fiber(&node->rightFiber, time);
	checkForDegrade_Fiber(&node->upFiber, time);
}

//************************************************//
//************************************************//

int findBindTime(Parameters* parameters, int time) {
	return (int)(time - log(urcw1_()) / (parameters->bindingRate * parameters->bindingSites) - parameters->timeStep / 2);
}

#ifdef NEWTIMELIST
void deleteOneMolecule(DynArray** timeList, int* totalMolecules, int* index) {
	(*totalMolecules)--;
	delAtIndexDA(timeList, index);
}
#else
void deleteOneMolecule(TimeList** timeList, int* totalMolecules, int* index) {
	(*totalMolecules)--;
	deleteOneTime(timeList, index);
}
#endif

void addUnboundToFiber(Fiber* fiber, Parameters* parameters, int time) {
	if (fiber->fiberIsDegraded) {
		fiber->unboundMolecules++;
	}
	else {
		fiber->unboundMolecules++;
		int bindTime = findBindTime(parameters, time);
#ifdef NEWTIMELIST
		addToDynArray(&fiber->unBoundTimeList, bindTime);
#else
		addToTimeList(&fiber->unBoundTimeList, bindTime);
#endif
	}
}

void addMultipleUnboundToFiber(Fiber* fiber, Parameters* parameters, int time, int count) {
	for (int iii = 0; iii < count; iii++) {
		addUnboundToFiber(fiber, parameters, time);
	}
}

//************************************************//
//************************************************//

void checkForUnbinds_Fiber(Fiber* fiber, Parameters* parameters, int time) {

	if (fiber->boundMolecules > 0) {
		int nextTime;
		int index = -1;
		int totalMolecules = fiber->boundMolecules;
		for (int iii = 0; iii < totalMolecules; iii++){
#ifdef NEWTIMELIST
			nextTime = retDynArray(fiber->boundTimeList, &index);
#else
			nextTime = retrieveNextTime(fiber->boundTimeList, &index);
#endif
			if (nextTime <= time) {
				deleteOneMolecule(&fiber->boundTimeList, &fiber->boundMolecules, &index);
				addUnboundToFiber(fiber, parameters, time);
			}
		}
	}
}

void checkForUnbinds_Node(Node* node, Parameters* parameters, int time) {
	checkForUnbinds_Fiber(&node->forwardFiber, parameters, time);
	checkForUnbinds_Fiber(&node->rightFiber, parameters, time);
	checkForUnbinds_Fiber(&node->upFiber, parameters, time);
}

//************************************************//
//************************************************//

bool moveMolecule(NodeGrid* grid, Parameters* parameters, int nodeIndex, int fiberDirection, int time) { //returns true if move is successful
	int newFiber = (int)(urcw1_() * 8);

	if (FORWARD == fiberDirection) {
		if (0 == newFiber || 2 == newFiber) { // move to up fiber or down fiber
			Fiber* fiber = &grid->nodeArray[nodeIndex].upFiber;
			addUnboundToFiber(fiber, parameters, time);
			return true;
		}
		else if (1 == newFiber) { //move to right fiber
			if (grid->rightBoundary && 0 == (nodeIndex + 1) % grid->nodeWidth) { //the farthest fiber on the right is not part of the simulation
				return false;
			}
			else{
				Fiber* fiber = &grid->nodeArray[nodeIndex].rightFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
		}
		else if (3 == newFiber) { //move to left fiber (aka right fiber on node directly to left
#if NEWCODE
				if (0 == nodeIndex % grid->nodeWidth) { //fiber on left is on next process to left, need to transfer
					if (grid->leftBoundary) {
						return false;
					}
					else {
						int bufferIndex = nodeIndex / grid->nodeWidth;
						grid->transferUnboundToLeftArray[bufferIndex]++;
						return true;
					}
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex - 1].rightFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
		}
#else
			if (grid->leftBoundary && 0 == nodeIndex % grid->nodeWidth) { //no fiber on the left
				return false;
			}
			else {
				if (0 == nodeIndex % grid->nodeWidth) { //fiber on left is on next process to left, need to transfer
					int bufferIndex = nodeIndex / grid->nodeWidth;
					grid->transferUnboundToLeftArray[bufferIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex - 1].rightFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
		}
#endif
		else if (4 == newFiber || 6 == newFiber) { //move to up fiber on node directly forward
			//don't need to worry about very top fiber since it is not part of simulation
			if (nodeIndex > (grid->totalNodesWithGhostRows - grid->nodeWidth - 1)) { //top row, need to transfer
				int bufferIndex = nodeIndex - (grid->totalNodesWithGhostRows - grid->nodeWidth);
				bufferIndex *= 2; //array is laid out in up, right, up, right, up, right...
				grid->transferUnboundToUpArray[bufferIndex]++;
				return true;
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex + grid->nodeWidth].upFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
		}
		else if (5 == newFiber) { //move to right fiber on node directly forward
			if (grid->rightBoundary && 0 == (nodeIndex + 1) % grid->nodeWidth){ //fiber on very far right is not part of simulation
				//|| grid->topBoundary && nodeIndex > grid->totalNodesWithGhostRows - grid->nodeWidth - 1) { //very top row, no subgrids above
				return false;
			}
			else {
				if (nodeIndex > (grid->totalNodesWithGhostRows - grid->nodeWidth - 1)) {
					int bufferIndex = nodeIndex - (grid->totalNodesWithGhostRows - grid->nodeWidth);
					bufferIndex *= 2; //array is laid out: up, right, up, right, up, right...
					bufferIndex++;
					grid->transferUnboundToUpArray[bufferIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex + grid->nodeWidth].rightFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
		}
		else if (7 == newFiber) { //move to left fiber on node directly forward (right fiber on node to up and left)
			if (grid->leftBoundary && 0 == nodeIndex % grid->nodeWidth) { //no fiber on the left
				return false;
			}
			else {
				if (nodeIndex == grid->totalNodesWithGhostRows - grid->nodeWidth) {
					grid->transferToUpperLeft++;
					return true;
				}
				else if (nodeIndex > (grid->totalNodesWithGhostRows - grid->nodeWidth - 1)) { //top row
					int bufferIndex = nodeIndex - (grid->totalNodesWithGhostRows - grid->nodeWidth) - 1;
					bufferIndex *= 2; //array is laid out: up, right, up, right, up, right...
					bufferIndex++;
					grid->transferUnboundToUpArray[bufferIndex]++;
					return true;
				}
				else if (0 == nodeIndex % grid->nodeWidth) { //fiber on left is on next process to left, need to transfer
					int bufferIndex = nodeIndex / grid->nodeWidth + 1;
					grid->transferUnboundToLeftArray[bufferIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex + grid->nodeWidth - 1].rightFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
		}
	}
	else if (UP == fiberDirection) { 
		if (0 == newFiber || 4 == newFiber) { //move to forward fiber on same node
			if (grid->topBoundary && nodeIndex > (grid->totalNodesWithGhostRows - grid->nodeWidth - 1)) { //very top fiber, not part of simulation
				return false;
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex].forwardFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
		}
		else if (1 == newFiber || 5 == newFiber) { //move to right fiber on same node
			if (grid->rightBoundary && 0 == (nodeIndex + 1) % grid->nodeWidth) { //very right fiber, not part of simulation
				return false;
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex].rightFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
		}
		else if (2 == newFiber || 6 == newFiber) { //move to back fiber (aka forward fiber of node directly back)
#if NEWCODE		
			if (nodeIndex < grid->nodeWidth) { //bottom row of process
				if (grid->bottomBoundary) {
					return false;
				}
				else {
					grid->transferUnboundToDownArray[nodeIndex]++;
					return true;
				}
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex - grid->nodeWidth].forwardFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
#else
			if (grid->bottomBoundary && nodeIndex < grid->nodeWidth) { //very last ghost row
				return false;
			}
			else {
				if (nodeIndex < grid->nodeWidth) { //bottom row of process
					grid->transferUnboundToDownArray[nodeIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex - grid->nodeWidth].forwardFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
#endif
		}
		else if (3 == newFiber || 7 == newFiber) { //move to left fiber (right fiber on node to the left)
#if NEWCODE
			if (0 == nodeIndex % grid->nodeWidth) { //on the left, need to transfer to process on left
				if (grid->leftBoundary) {
					return false;
				}
				else {
					int bufferIndex = nodeIndex / grid->nodeWidth;
					grid->transferUnboundToLeftArray[bufferIndex]++;
					return true;
				}
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex - 1].rightFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
#else
			if (grid->leftBoundary && 0 == nodeIndex % grid->nodeWidth) {
				return false;
			}
			else {
				if (0 == nodeIndex % grid->nodeWidth) { //on the left, need to transfer to process on left
					int bufferIndex = nodeIndex / grid->nodeWidth;
					grid->transferUnboundToLeftArray[bufferIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex - 1].rightFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
#endif
		}
	}
	else if (RIGHT == fiberDirection) {
		if (0 == newFiber || 2 == newFiber) { //move to up (0) or down (2) fiber to the left (on same node)
			Fiber* fiber = &grid->nodeArray[nodeIndex].upFiber;
			addUnboundToFiber(fiber, parameters, time);
			return true;
		}
		else if (1 == newFiber) { //move to back fiber (forward fiber of node directly back)
#if NEWCODE
			if (nodeIndex < grid->nodeWidth) {
				if (grid->bottomBoundary) {
					return false;
				}
				else {
					//int bufferIndex = nodeIndex - (grid->totalNodesWithGhostRows - grid->nodeWidth);
					int bufferIndex = nodeIndex;
					grid->transferUnboundToDownArray[bufferIndex]++;
					return true;
				}
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex - grid->nodeWidth].forwardFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
#else
			if (grid->bottomBoundary && nodeIndex < grid->nodeWidth) { //very bottom ghost row
				return false;
			}
			else {
				if (nodeIndex < grid->nodeWidth) {
					//int bufferIndex = nodeIndex - (grid->totalNodesWithGhostRows - grid->nodeWidth);
					int bufferIndex = nodeIndex;
					grid->transferUnboundToDownArray[bufferIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex - grid->nodeWidth].forwardFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
#endif
		}
		else if (3 == newFiber) { //move to forward fiber
			if (grid->topBoundary && nodeIndex > (grid->totalNodesWithGhostRows - grid->nodeWidth - 1)) { // the forward fiber on the very top row is not part of the simulation
				return false;
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex].forwardFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
		}
		else if (4 == newFiber || 6 == newFiber) { //move to up fiber on node to directly to right
			//don't need to worry about fiber all the way to right since it's not part of the simulation
			if (0 == (nodeIndex + 1) % grid->nodeWidth) { //if the node is on the right of the process
				int bufferIndex = nodeIndex / grid->nodeWidth;
				bufferIndex *= 2;
				bufferIndex++;
				grid->transferUnboundToRightArray[bufferIndex]++;
				return true;
			}
			else {
				Fiber* fiber = &grid->nodeArray[nodeIndex + 1].upFiber;
				addUnboundToFiber(fiber, parameters, time);
				return true;
			}
		}
		//add code
		//need to add code to handle move to right when on right fiber on right side of sub grid
		else if (5 == newFiber) { //move to forward fiber on node to down right
#if NEWCODE
			if (nodeIndex < grid->nodeWidth) { //bottom row
				if (grid->bottomBoundary) {
					return false;
				}
				else if (grid->nodeWidth == nodeIndex + 1) {
					grid->transferToLowerRight++;
					return true;
				}
				else {
					int bufferIndex = nodeIndex + 1;
					grid->transferUnboundToDownArray[bufferIndex]++;
					return true;
				}
			}
			else { //rows above bottom row
				if (0 == (nodeIndex + 1) % grid->nodeWidth) { //if on the right of the subgrid, move to subgrid to the right
					int bufferIndex = nodeIndex / grid->nodeWidth - 1;
					bufferIndex *= 2;
					grid->transferUnboundToRightArray[bufferIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex - grid->nodeWidth + 1].forwardFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
#else
			if (grid->bottomBoundary && nodeIndex < grid->nodeWidth) { //very bottom ghost row
				return false;
			}
			else {
				if (nodeIndex < grid->nodeWidth) { //bottom row
					if (grid->nodeWidth == nodeIndex + 1) {
						grid->transferToLowerRight++;
						return true;
					}
					else {
						int bufferIndex = nodeIndex + 1;
						grid->transferUnboundToDownArray[bufferIndex]++;
						return true;
					}
				}
				else { //rows above bottom row
					if (0 == (nodeIndex + 1) % grid->nodeWidth) { //if on the right of the subgrid, move to subgrid to the right
						int bufferIndex = nodeIndex / grid->nodeWidth - 1;
						bufferIndex *= 2;
						grid->transferUnboundToRightArray[bufferIndex]++;
						return true;
					}
					else {
						Fiber* fiber = &grid->nodeArray[nodeIndex - grid->nodeWidth + 1].forwardFiber;
						addUnboundToFiber(fiber, parameters, time);
						return true;
					}
				}
			}
#endif
		}
		else if (7 == newFiber) { //move to forward fiber on node directly to right
			//don't need to about the fiber farthest to the right since it's not part of the simulation
			if (grid->topBoundary && nodeIndex > (grid->totalNodesWithGhostRows - grid->nodeWidth - 1)) {
				return false;
			}
			else {
				if (0 == (nodeIndex + 1) % grid->nodeWidth) {
					int bufferIndex = nodeIndex / grid->nodeWidth;
					bufferIndex *= 2;
					grid->transferUnboundToRightArray[bufferIndex]++;
					return true;
				}
				else {
					Fiber* fiber = &grid->nodeArray[nodeIndex + 1].forwardFiber;
					addUnboundToFiber(fiber, parameters, time);
					return true;
				}
			}
		}
	}
	printf("File: %s Line: %d Should not reach bottom of moveMolecule()! Aborting\n", __FILE__ , __LINE__);
	fflush(stdout);
	MPI_Abort(MPI_COMM_WORLD, 100);
}

//************************************************//
//************************************************//

void setDegradeTime(Fiber* fiber, InputData* inputData, Parameters* parameters, int time, int colr2) {
	double random = urcw1_();
	int r400 = ceil(random * 100) + 1; //needs to be changed, probably
	if (r400 <= inputData->lysesPerBlock[colr2 - 1]) {
		double slope = inputData->lysisTime[r400][colr2 - 1] - inputData->lysisTime[r400 - 1][colr2 - 1];
		double step = (random * 100) - (r400 - 1); //needs to be changed, probably
		double timeToDegrade = inputData->lysisTime[r400][colr2 - 1] + step * slope;
		int newTime = (int)(time + floor(timeToDegrade / parameters->timeStep));
		if (newTime < fiber->degradeTime) {
			fiber->degradeTime = newTime;
		}
	}
}

int findUnbindTime(Fiber* fiber, InputData* inputData, Parameters* parameters, int time) {
	double random = urcw1_();
	int colr2 = ceil(random * 100); //needs to be changed, probably
	double slope = inputData->unbindingTimeDistribution[colr2] - inputData->unbindingTimeDistribution[colr2 - 1];
	double step = (random * 100) - colr2; //needs to be changed, probably
	double timeToUnbind = inputData->unbindingTimeDistribution[colr2 - 1] + (step*slope);
	int newTime = time + timeToUnbind / parameters->timeStep;
#ifdef NEWTIMELIST
	addToDynArray(&fiber->boundTimeList, newTime);
#else
	addToTimeList(&fiber->boundTimeList, newTime);
#endif
	return colr2;
}

void bindMolecule(Fiber* fiber, InputData* inputData, Parameters* parameters, int time) {
	fiber->boundMolecules++;
	int colr2 = findUnbindTime(fiber, inputData, parameters, time);
	setDegradeTime(fiber, inputData, parameters, time, colr2);
}

//************************************************//
//************************************************//

void unboundPath_Fiber(NodeGrid* grid, int nodeIndex, Fiber* fiber, int fiberDirection, Parameters* parameters, InputData* inputData, int time) {

	if (fiber->unboundMolecules > 0) {
		if (fiber->fiberIsDegraded) { //special class of unbound molecules on a degraded fiber, can only stay or move
			//when fiber is degraded, all molecules are unbound and there is no timelist
			int total = fiber->unboundMolecules;
			for (int iii = 0; iii < total; iii++) {
				double random = urcw1_();
				//if (!(random <= 1 - parameters->movingProbability)) { //may not need to invert, but it more closely resembles Dr. Paynter's code
				if (random <= parameters->movingProbability) {
					if (moveMolecule(grid, parameters, nodeIndex, fiberDirection, time)) {
						fiber->unboundMolecules--;
					}
				}
			}
		}
		else { //fiber is not degraded
			int index = -1;
			int total = fiber->unboundMolecules; //unBound molecules will be changing, so we need to use a different variable for the loop

			for (int iii = 0; iii < total; iii++) {
#ifdef NEWTIMELIST
				int bindTime = retDynArray(fiber->unBoundTimeList, &index);
#else
				int bindTime = retrieveNextTime(fiber->unBoundTimeList, &index);
#endif
				double random = urcw1_();
				if (random < 1 - parameters->movingProbability) {
					//molecule will bind or stay unbound at same fiber
					if (bindTime <= time) {
						bindMolecule(fiber, inputData, parameters, time);
						deleteOneMolecule(&fiber->unBoundTimeList, &fiber->unboundMolecules, &index);
					}
				}
				else {
					if (bindTime <= time) {
						if (urcw1_() > (time - bindTime / parameters->timeStep)) {
							if (moveMolecule(grid, parameters, nodeIndex, fiberDirection, time)) {
								deleteOneMolecule(&fiber->unBoundTimeList, &fiber->unboundMolecules, &index);
							}
						}
						else {
							bindMolecule(fiber, inputData, parameters, time);
							deleteOneMolecule(&fiber->unBoundTimeList, &fiber->unboundMolecules, &index);
						}
					}
					else { //bindTime is greater than current time
						if (moveMolecule(grid, parameters, nodeIndex, fiberDirection, time)) {
							deleteOneMolecule(&fiber->unBoundTimeList, &fiber->unboundMolecules, &index);
						}
					}
				}
			} //end loop
		}
	}
}

void unboundPath_Node(NodeGrid* grid, int nodeIndex, Node* node, Parameters* parameters, InputData* inputData, int time) {
	unboundPath_Fiber(grid, nodeIndex, &node->forwardFiber, FORWARD, parameters, inputData, time);
	unboundPath_Fiber(grid, nodeIndex, &node->rightFiber, RIGHT, parameters, inputData, time);
	unboundPath_Fiber(grid, nodeIndex, &node->upFiber, UP, parameters, inputData, time);
}

//************************************************//
//************************************************//

void simulationLoop(NodeGrid* grid, Parameters* parameters, InputData* inputData, RankGrid* rankGrid, Files* files) {

	int captureInterval = parameters->totalTimeSteps / (parameters->totalTime / 10);
	int nextCapture = captureInterval;

	for (int time = 0; time < parameters->totalTimeSteps; time++) {
		for (int nodeIndex = 0; nodeIndex < grid->totalNodesWithGhostRows; nodeIndex++) {
			if (nodeIndex > grid->indexOfLastGhostNode) {
				checkForDegrade_Node(&grid->nodeArray[nodeIndex], time);
				checkForUnbinds_Node(&grid->nodeArray[nodeIndex], parameters, time);
			}
			unboundPath_Node(grid, nodeIndex, &grid->nodeArray[nodeIndex], parameters, inputData, time);
		}
		transferToOtherRanks(grid, rankGrid, parameters, time);
		if (time == nextCapture) {
			printBoundMoleculesToFile(grid, rankGrid, files->boundMolecules);
			printUnboundMoleculesToFile(grid, rankGrid, files->unboundMolecules);
			printIsDegradedToFile(grid, rankGrid, files->fibers);
			nextCapture += captureInterval;
		}
	}
	debug_checkAllMolecules(grid, parameters, __FILE__, __LINE__);
}

//************************************************//
//************************************************//

int debug_countLocalMolecules(NodeGrid* grid) {
	int count = 0;

	for (int iii = 0; iii < grid->totalNodesWithGhostRows; iii++) {
		count += (grid->nodeArray[iii].forwardFiber.unboundMolecules + grid->nodeArray[iii].forwardFiber.boundMolecules);
		count += (grid->nodeArray[iii].upFiber.unboundMolecules + grid->nodeArray[iii].upFiber.boundMolecules);
		count += (grid->nodeArray[iii].rightFiber.unboundMolecules + grid->nodeArray[iii].rightFiber.boundMolecules);
	}

	if (NULL != grid->transferUnboundToUpArray) {
		for (int iii = 0; iii < grid->transferUnboundToUpArrayCount; iii++) {
			count += grid->transferUnboundToUpArray[iii];
		}
	}
	if (NULL != grid->transferUnboundToRightArray) {
		for (int iii = 0; iii < grid->transferUnboundToRightArrayCount; iii++) {
			count += grid->transferUnboundToRightArray[iii];
		}
	}
	if (NULL != grid->transferUnboundToDownArray) {
		for (int iii = 0; iii < grid->transferUnboundToDownArrayCount; iii++) {
			count += grid->transferUnboundToDownArray[iii];
		}
	}
	if (NULL != grid->transferUnboundToLeftArray) {
		for (int iii = 0; iii < grid->transferUnboundToLeftArrayCount; iii++) {
			count += grid->transferUnboundToLeftArray[iii];
		}
	}

	return count;
}

void debug_checkAllMolecules(NodeGrid* grid, Parameters* parameters, char* file, int line) {

	int count = debug_countLocalMolecules(grid);

	int total;
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&count, &total, 1, MPI_INT, MPI_SUM, HEAD_RANK, MPI_COMM_WORLD);
	if (HEAD_RANK == grid->MPIRank && total != parameters->totalMolecules) {
		printf("Molecule count: %i does not match in file: %s line:%i\n", total, file, line);
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, 0);
	}
}





