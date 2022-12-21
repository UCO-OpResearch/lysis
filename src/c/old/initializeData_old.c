#include "all.h"

bool loadLysisTimeFromFile(Files* files, InputData* inputData, Parameters* parameters) {

	FILE* file;

	file = fopen(files->lysisTimeFile, "r");
	printf("0\n");
	if (NULL == file) {
		printf("intializeData.c: loadLysisTimeFromFile: lysis time file %s failed to open.\n", files->lysisTimeFile);
		return false;
	}
	else {
		double** lysisTime = inputData->lysisTime;

		char* line = NULL;
		size_t length = 0;
		printf("1\n");
		int iii = 0;
		while (EOF != getline(&line, &length, file)){
			if (iii == 1) { printf("first loop"); }
			int jjj = 0;
			while (EOF != sscanf(line, "%f", &lysisTime[iii][jjj])) {
				jjj++;
				if (jjj == 1) { printf("second loop"); }
			}

			if (jjj != parameters->lysisBlocks) {
				printf("initializeData.c: loadLysisTimeFromFile: expected columns (lysis blocks) is %i. Number found is %i on row %i of file %s.\n", parameters->lysisBlocks, jjj, iii, files->lysisTimeFile);
				return false;
			}

			iii++;
		}
		printf("2\n");
		if (iii != parameters->maxLysesPerBlock) {
			printf("intializeData.c: loadLysisTimeFromFile: expected rows (maxLysesPerBlock) is %i. Found rows is %i in file  %s.\n", parameters->maxLysesPerBlock, iii, files->lysisTimeFile);
			return false;
		}

		if (0 == feof(file) || 0 != ferror(file)) {
			printf("initializeData.c: loadLysisTimeFromFile: error loading data\n");
			return false;
		}
		printf("3\n");
		free(line);
		fclose(file);
		printf("4\n");
		return true;
	}
}

bool loadUnbindingTimeFromFile(Files* files, InputData* inputData, Parameters* parameters) {

	FILE* file;
	file = fopen(files->unbindingTimeFile, "r");
	
	double* unbindingTimeDistribution = inputData->UnbindingTimeDistribution;

	if (NULL == file) {
		printf("initializeData.c: loadUnbindingTimeFromFile: unbinding time file %s failed to open.\n", files->unbindingTimeFile);
		return false;
	}
	else {

		int iii = 0;
		while (EOF != fscanf(file, "%f", &unbindingTimeDistribution[iii])) {
			iii++;
		}

		if (iii != parameters->unbindingTimeDist) {
			printf("initializeData.c: loadUnbindingTimeFromFile: entries mismatch, expected entries = %i, loaded entries = %i", parameters->unbindingTimeDist, iii);
			return false;
		}
	}

	fclose(file);
	return true;
}

bool loadLysesPerBlockFromFile(Files* files, InputData* inputData, Parameters* parameters) {

	FILE* file;
	file = fopen(files->totalLysesFile, "r");

	int* lysesPerBlock = inputData->lysesPerBlock;

	if (NULL == file) {
		printf("initializeData.c: loadLysesPerBlockFromFile: total lyses file %s failed to open.\n", files->totalLysesFile);
		return false;
	}
	else {

		int iii = 0;
		while (EOF != fscanf(file, "%i", &lysesPerBlock[iii])) {
			iii++;
		}

		if (iii != parameters->lysisBlocks) {
			printf("initializeData.c: loadLysesPerBlockFromFile: entries mismatch, expected entries = %i, loaded entries = %i", parameters->lysisBlocks, iii);
			return false;
		}
	}

	fclose(file);
	return true;
}

bool loadInputData(Files* files, InputData* inputData, Parameters* parameters) {

	if (!loadLysisTimeFromFile(files, inputData, parameters)) {
		return false;
	}
	/*
	if (!loadUnbindingTimeFromFile(files, inputData, parameters)) {
		return false;
	}

	if (!loadLysesPerBlockFromFile(files, inputData, parameters)) {
		return false;
	}
	*/
	return true;
}