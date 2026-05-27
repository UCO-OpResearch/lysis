#include "all.h"

#define ELEMENTS 40
#define ARRAYS ELEMENTS
#define REDUCE 10 // sets at what point the arrays are reduced

void addToDynArray(DynArray** in, int newTime) {

	DynArray* d = *in;

	if (NULL == d) {
		d = (DynArray*)malloc(sizeof(DynArray));
		d->total = 1;
		d->potentials = ARRAYS;
		d->allocations = 1;

		d->counts = (int*)malloc(sizeof(int) * ARRAYS);
		d->counts[0] = 1;
		for (int iii = 1; iii < ARRAYS; iii++) {
			d->counts[iii] = 0;
		}

		d->arrays = (int**)malloc(sizeof(int*) * ARRAYS);
		d->arrays[0] = (int*)malloc(sizeof(int) * ELEMENTS);
		for (int iii = 1; iii < ARRAYS; iii++) {
			d->arrays[iii] = NULL;
		}
		d->arrays[0][0] = newTime;
		for (int iii = 1; iii < ELEMENTS; iii++) {
			d->arrays[0][iii] = -1;
		}
		*in = d;
	}
	else {
		if (d->total == d->allocations * ELEMENTS) { // if true, arrays are full, need to add more
			int index;
			if (d->allocations == d->potentials) { // no free pointers-to-arrays need to expand
				d->potentials = d->potentials + ARRAYS;
				int** toDelete = d->arrays;
				int* oldCounts = d->counts;
				d->arrays = (int**)malloc(sizeof(int*) * d->potentials);
				d->counts = (int*)malloc(sizeof(int) * d->potentials);
				for (int iii = 0; iii < d->allocations; iii++) {
					d->arrays[iii] = toDelete[iii];
					d->counts[iii] = ELEMENTS;
				}
				for (int iii = d->allocations; iii < d->potentials; iii++) {
					d->arrays[iii] = NULL;
					d->counts[iii] = 0;
				}
				index = d->allocations;
				free(toDelete);
				free(oldCounts);
			}
			else {// otherwise just find the null ptr
				index = 0;
				while (d->arrays[index] != NULL) {
					index++;
				}
			}

			d->allocations++;
			d->arrays[index] = (int*)malloc(sizeof(int) * ELEMENTS);
			for (int iii = 0; iii < ELEMENTS; iii++) {
				d->arrays[index][iii] = -1;
			}

			d->total++;
			d->counts[index]++;
			d->arrays[index][0] = newTime;
		}
		else { // arrays are not full, find free element and insert newTime
			int index = 0;
			while (d->counts[index] == ELEMENTS || NULL == d->arrays[index]) {
				index++;
			}
			int iii = 0;
			while (-1 != d->arrays[index][iii]) {
				iii++;
			}
			d->total++;
			d->counts[index]++;
			d->arrays[index][iii] = newTime;
		}
	}
}

int retDynArray(DynArray* d, int* index) {

	int temp = *index;

	int which = temp / ELEMENTS;
	int element = temp % ELEMENTS;

	element++;
	temp++;
	if (ELEMENTS == element) { //if true, need to go to next array
		which++;
		while (NULL == d->arrays[which]) { // skip over any null ptrs
			which++;
			temp += ELEMENTS; // null ptr is considered to be an empty array
		}
		element = 0;
	}

	while (0 == d->counts[which]) {
		temp += ELEMENTS;
		which++;
	}

	while (-1 == d->arrays[which][element]) {
		element++;
		temp++;
		if (ELEMENTS == element) {
			which++;
			while (NULL == d->arrays[which]) { // skip over any null ptrs
				which++;
				temp += ELEMENTS;
			}
			element = 0;
		}
	}

	*index = temp;
	return d->arrays[which][element];

}

void delAtIndexDA(DynArray** in, int* inIndex) {

	int index = *inIndex;
	DynArray* d = *in;

	int which = index / ELEMENTS;
	int element = index % ELEMENTS;

	if (-1 == d->arrays[which][element]) {
		fprintf(stderr, "Logical error in delAtIndexDA function in %s. Trying to delete already deleted value.",  __FILE__);
	}

	d->arrays[which][element] = -1;
	d->counts[which]--;
	d->total--;

	if (0 == d->counts[which]) {
		free(d->arrays[which]);
		d->arrays[which] = NULL;
		d->allocations--;
		index += ELEMENTS - 1 - element; //move to end of array
		if (0 == d->allocations) {
			free(d->arrays);
			free(d);
			*in = NULL;
		}
		else if (d->potentials - d->allocations == REDUCE + ARRAYS) {
			int newP = d->potentials - ARRAYS;
			int** newArrays = (int**)malloc(sizeof(int*) * newP);
			int nnn = 0;
			for (int iii = 0; iii < d->potentials; iii++) {
				if (NULL != d->arrays[iii]) {
					newArrays[nnn] = d->arrays[iii];
					nnn++;
				}
				else {
					if (iii < which + 1) {
						index -= ELEMENTS;
					}
				}
			}
			while (nnn < newP) {
				newArrays[nnn] = NULL;
				nnn++;
			}
			d->potentials = newP;
			free(d->arrays);
			d->arrays = newArrays;
		}
	}
	*inIndex = index;
}

void delEntireDA(DynArray** in) {
	DynArray* d = *in;
	if (NULL != d->arrays) {
		for (int iii = 0; iii < d->potentials; iii++) {
			free(d->arrays[iii]);
		}
		free(d->arrays);
		free(d);
		*in = NULL;
	}

}
