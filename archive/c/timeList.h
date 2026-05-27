#pragma once

#define TIME_BLOCK_SIZE 48 //sets how many ints are allocated at a time for time lists. Also used in timeList functions

typedef struct TimeList TimeList;

struct TimeList {
	long times[TIME_BLOCK_SIZE];
	TimeList* next;
};

void addToTimeList(TimeList** timeList, int newTime);
void deleteOneTime(TimeList** timeList, int* index);
void deleteTimeList(TimeList** timeList);
int retrieveNextTime(TimeList* timeList, int* index);
