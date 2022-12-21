#include "all.h"

bool loadLysisTimeFromFile(Files* files, InputData* inputData, Parameters* parameters) {

	FILE* file;
	file = fopen(files->lysisTimeFile, "r");

	if (NULL == file) {
		printf("intializeData.c: loadLysisTimeFromFile: lysis time file %s failed to open.\n", files->lysisTimeFile);
		return false;
	}
	else {
		double** lysisTime = inputData->lysisTime;

		char* line = NULL;
		size_t length = 0;
		int iii = 0;
		while (EOF != getline(&line, &length, file)){

			int jjj = 0;
			FILE* ftemp = fmemopen(line, strlen(line), "r");

			while (EOF != fscanf(ftemp, "%lf", &inputData->lysisTime[iii][jjj])) {
				jjj++;
			}

			iii++;

			if (jjj != parameters->lysisBlocks) {
				printf("initializeData.c: loadLysisTimeFromFile: expected columns (lysis blocks) is %i. Number found is %i on row %i of file %s.\n", parameters->lysisBlocks, jjj, iii, files->lysisTimeFile);
				return false;
			}
		}

		if (iii != parameters->maxLysesPerBlock) {
			printf("intializeData.c: loadLysisTimeFromFile: expected rows (maxLysesPerBlock) is %i. Found rows is %i in file  %s.\n", parameters->maxLysesPerBlock, iii, files->lysisTimeFile);
			return false;
		}

		if (0 == feof(file) || 0 != ferror(file)) {
			printf("initializeData.c: loadLysisTimeFromFile: error loading data\n");
			return false;
		}

		free(line);
		fclose(file);

		//output2dArrayToFile(inputData, parameters);
		return true;
	}
}

bool loadUnbindingTimeFromFile(Files* files, InputData* inputData, Parameters* parameters) {

	FILE* file;
	file = fopen(files->unbindingTimeFile, "r");

	if (NULL == file) {
		printf("initializeData.c: loadUnbindingTimeFromFile: unbinding time file %s failed to open.\n", files->unbindingTimeFile);
		return false;
	}
	else {

		int iii = 0;
		double temp;
		while (EOF != fscanf(file, "%lf", &inputData->unbindingTimeDistribution[iii])) {
			iii++;
		}

		if (iii != parameters->unbindingTimeDist) {
			printf("initializeData.c: loadUnbindingTimeFromFile: entries mismatch, expected entries = %i, loaded entries = %i", parameters->unbindingTimeDist, iii);
			return false;
		}
	}

	//outputArrayToFile(inputData, parameters);

	fclose(file);
	return true;
}

bool loadLysesPerBlockFromFile(Files* files, InputData* inputData, Parameters* parameters) {

	FILE* file;
	file = fopen(files->totalLysesFile, "r");

	if (NULL == file) {
		printf("initializeData.c: loadLysesPerBlockFromFile: total lyses file %s failed to open.\n", files->totalLysesFile);
		return false;
	}
	else {

		int iii = 0;
		while (EOF != fscanf(file, "%i", &inputData->lysesPerBlock[iii])) {
			iii++;
		}

		if (iii != parameters->lysisBlocks) {
			printf("initializeData.c: loadLysesPerBlockFromFile: entries mismatch, expected entries = %i, loaded entries = %i", parameters->lysisBlocks, iii);
			return false;
		}
	}

	//outputLysesPerBlockToFile(inputData, parameters);

	fclose(file);
	return true;
}

bool loadInputData(Files* files, InputData* inputData, Parameters* parameters) {

	if (!loadLysisTimeFromFile(files, inputData, parameters)) {
		return false;
	}

	if (!loadUnbindingTimeFromFile(files, inputData, parameters)) {
		return false;
	}

	if (!loadLysesPerBlockFromFile(files, inputData, parameters)) {
		return false;
	}

	return true;
}

//for debugging
void output2dArrayToFile(InputData* inputData, Parameters* parameters) {
	FILE* output = fopen("output/test.dat", "w");

	for (int iii = 0; iii < parameters->maxLysesPerBlock; iii++) {
		for (int jjj = 0; jjj < parameters->lysisBlocks; jjj++) {
			fprintf(output, "   %.16e", inputData->lysisTime[iii][jjj]);
		}
		fprintf(output, "%s", "\r\n");
	}

	fclose(output);
}

void outputArrayToFile(InputData* inputData, Parameters* parameters) {
	FILE* output = fopen("output/test-tsect.dat", "w");

	for (int iii = 0; iii < parameters->unbindingTimeDist; iii++) {
		fprintf(output, "   %.16e\r\n", inputData->unbindingTimeDistribution[iii]);
	}

	fclose(output);
}

void outputLysesPerBlockToFile(InputData* inputData, Parameters* parameters) {
	FILE* output = fopen("output/test-lenlysisvect.dat", "w");

	for (int iii = 0; iii < parameters->lysisBlocks; iii++) {
		fprintf(output, "%i\r\n", inputData->lysesPerBlock[iii]);
	}

	fclose(output);
}