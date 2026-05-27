#pragma once

typedef struct DynArray DynArray;

struct DynArray {

	int total;
	int potentials;
	int allocations;
	int* counts;
	int** arrays;

};

void addToDynArray(DynArray** in, int newTime);
int retDynArray(DynArray* d, int* index);
void delAtIndexDA(DynArray** in, int* inIndex);
void delEntireDA(DynArray** in);