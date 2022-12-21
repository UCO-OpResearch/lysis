#include "all.h"

/*if (0 == nodeIndex % grid->nodeWidth) {
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
}*/

bool nextTimeList(TimeList** workingTL) {
	if (NULL != (*workingTL)->next) {
		(*workingTL) = (*workingTL)->next;
		return true;
	}
	else {
		return false;
	}
}

void test_printTimeList(TimeList** timeList) {
	TimeList* workingTL = (*timeList);
	if (NULL != workingTL) {
		int count = 0;
		do {
			printf("timeList: %i\n", count);
			for (int iii = 0; iii < TIME_BLOCK_SIZE; iii++) {
				printf("%i\n", workingTL->times[iii]);
			}
			count++;
		} while (nextTimeList(&workingTL));
	}
}

void test_printTimeListToFile(TimeList** timeList, char* fileName) {
	FILE* file = fopen(fileName, "a");
	TimeList* workingTL = (*timeList);

	if (NULL != workingTL) {
		int count = 0;
		do {
			for (int iii = 0; iii < TIME_BLOCK_SIZE; iii++) {
				if(-1 != workingTL->times[iii])
				fprintf(file, "%i,", workingTL->times[iii]);
			}
		} while (nextTimeList(&workingTL));
		fprintf(file, "\r\n");
	}
	fclose(file);
}

/*void test_memoryOperations() {
	TimeList* timeList = NULL;
	//printf("%i\n", (int)sizeof(timeList));
	//printf("%i\n", (int)sizeof(long));

	int timeBlocks = 5;
	int totalTimes = timeBlocks * TIME_BLOCK_SIZE;
	for (int iii = 0; iii < totalTimes; iii++) {
		addToTimeList(&timeList, iii);
	}

	int index;
	for (int iii = 16; iii < 24; iii++) { //delete middle
		index = iii;
		deleteOneTime(&timeList, &totalTimes, &index);
	}

	test_printTimeList(&timeList);

	printf("index: %i\n", index);
	*/

	/****************************************/

	/*TimeList* timeList2 = NULL;

	for (int iii = 0; iii < totalTimes; iii++) {
		addToTimeList(&timeList2, iii);
	}

	for (int iii = 32; iii < 40; iii++) { //delete last
		index = iii;
		deleteOneTime(&timeList2, &totalTimes, &index);
	}

	//test_printTimeList(&timeList2);

	//printf("index: %i\n", index);*/

	/****************************************/

	/*TimeList* timeList3 = NULL;

	for (int iii = 0; iii < totalTimes; iii++) {
		addToTimeList(&timeList3, iii);
	}

	for (int iii = 0; iii < 8; iii++) { //delete first
		index = iii;
		deleteOneTime(&timeList3, &totalTimes, &index);
	}

	//test_printTimeList(&timeList3);

	//printf("index: %i\n", index);*/

	/****************************************/

	/*TimeList* timeList4 = NULL;

	for (int iii = 0; iii < totalTimes; iii++) {
		addToTimeList(&timeList4, iii);
	}

	//test_printTimeList(&timeList4);

	index = -1;
	while(index < totalTimes - 1) {
		printf("retrieve: %i\n", retrieveNextTime(timeList4, &index));
	}

	deleteTimeList(&timeList4);*/

	/****************************************/

	/*int timeBlocks = 5;
	int totalTimes = timeBlocks * TIME_BLOCK_SIZE;

	TimeList* timeList5 = NULL;

	for (int iii = 0; iii < totalTimes; iii++) {
		addToTimeList(&timeList5, iii);
	}

	int deletions = 0.50 * totalTimes;
	for (int iii = 0; iii < deletions; iii++) {
		int toDelete = (int)(urcw1_() * totalTimes);
		deleteOneTime(&timeList5, &totalTimes, &toDelete);
	}

	test_printTimeList(&timeList5);

	int index = -1;
	while (index < totalTimes - 1) {
		printf("retrieve: %i\n", retrieveNextTime(timeList5, &index));
	}

	int loopCount = 0;
	while (NULL != timeList5) {
		int memBlocks = 1;
		TimeList* workingTL = timeList5;
		while (NULL != workingTL->next) {
			memBlocks++;
			workingTL = workingTL->next;
		}
		int toDelete = (int)(urcw1_() * memBlocks * TIME_BLOCK_SIZE);
		deleteOneTime(&timeList5, &totalTimes, &toDelete);
		loopCount++;
	}*/

	/****************************************/

	/*int timeBlocks = 5;
	int totalTimes = timeBlocks * TIME_BLOCK_SIZE;

	TimeList* timeList = NULL;

	for (int iii = 0; iii < totalTimes; iii++) {
		addToTimeList(&timeList, iii);
	}

	int loopCount = 0;
	while (NULL != timeList && loopCount < 500000) {
		if (urcw1_() > 0.50) {
			int memBlocks = 1;
			TimeList* workingTL = timeList;
			while (NULL != workingTL->next) {
				memBlocks++;
				workingTL = workingTL->next;
			}
			int toDelete = (int)(urcw1_() * memBlocks * TIME_BLOCK_SIZE);
			deleteOneTime(&timeList, &toDelete);
		}
		else {
			addToTimeList(&timeList, loopCount);
		}
		loopCount++;
	}
}*/