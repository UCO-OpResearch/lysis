#include "all.h"

int retrieveNextTime(TimeList* timeList, int* index) { //retrieves next valid time after given index

	int workingInd = *index + 1;

	int structIndex = workingInd / TIME_BLOCK_SIZE;
	int timeArrayIndex = workingInd % TIME_BLOCK_SIZE;

	TimeList* workingTL = timeList;
	for (int iii = 0; iii < structIndex; iii++) {
		workingTL = workingTL->next;
	}

	while (-1 == workingTL->times[timeArrayIndex]) {
		if (TIME_BLOCK_SIZE - 1 == timeArrayIndex) {
			workingTL = workingTL->next; //shouldn't have to check for null ptr since retrieveNextTime() shouldn't be called when there isn't a next time value
			timeArrayIndex = 0;
		}
		else {
			timeArrayIndex++;
			workingInd++;
		}
	}

	*index = workingInd;
	return 	workingTL->times[timeArrayIndex];
}

void deleteOneTime(TimeList** timeList, int* index) {// deletes time at index

	int structIndex = *index / TIME_BLOCK_SIZE;
	int timeArrayIndex = *index % TIME_BLOCK_SIZE;

	TimeList* workingTL = (*timeList);
	TimeList* previousTL = NULL;
	for (int iii = 0; iii < structIndex; iii++) {
		previousTL = workingTL;
		workingTL = workingTL->next;
	}
	workingTL->times[timeArrayIndex] = -1;

	bool occupied = false;
	for (int iii = 0; iii < TIME_BLOCK_SIZE && !occupied; iii++) {
		if (-1 != workingTL->times[iii]) {
			occupied = true;
		}
	}

	if (!occupied) {
		if (0 == structIndex) {
			if (NULL == workingTL->next) { // first and only timeList
				free(workingTL);
				*timeList = NULL;
			}
			else { // workingTL is first but not last
				*timeList = workingTL->next;
				free(workingTL);
			}
		}
		else {
			if (NULL == workingTL->next) { // end of list
				previousTL->next = NULL;
				free(workingTL);
			}
			else { // timeList is in the middle
				previousTL->next = workingTL->next;
				free(workingTL);
			}
		}
		*index = structIndex * TIME_BLOCK_SIZE - 1;
	}
}

void deleteTimeList(TimeList** timeList) {
	TimeList* next;
	while (NULL != (*timeList)->next) {
		next = (*timeList)->next;
		free((*timeList));
		(*timeList) = next;
	}
	free(*timeList);
	(*timeList) = NULL;
}

void initializeTimeListStruct(TimeList** timeList, int newTime) {
	(*timeList) = (TimeList*)malloc(sizeof(TimeList));
	(*timeList)->times[0] = newTime;
	for (int iii = 1; iii < TIME_BLOCK_SIZE; iii++) {
		(*timeList)->times[iii] = -1;
	}
	(*timeList)->next = NULL;
}

void addToTimeList(TimeList** timeList, int newTime) { //pass in old moleculeCount, function assumes count hasn't been incremented yet
	if (NULL == (*timeList)) {
		initializeTimeListStruct(timeList, newTime);
	}
	else {
		TimeList* workingTL = (*timeList);
		while (NULL != workingTL->next) {
			workingTL = workingTL->next;
		}
		if (-1 == workingTL->times[TIME_BLOCK_SIZE - 1]) { //timeList is not full
			int index = TIME_BLOCK_SIZE - 2;
			while (-1 == workingTL->times[index]) { //index should not get to -1 since an empty times array should not exist
				index--;
			}
			workingTL->times[index + 1] = newTime;
		}
		else { //timelist is full
			initializeTimeListStruct(&workingTL->next, newTime);
		}
	}
}