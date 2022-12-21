#pragma once

#define HEAD_RANK 0
#define FIBERS_PER_NODE 3

typedef struct {

	double bindingRate;
	double poreSize;
	double diffusionCoeff;
	double bindingSites; 
	double gridSymmetryDistance;
	int totalMolecules;
	double movingProbability;
	int totalTime;
	double timeStep;
	int totalTimeSteps;

	int unbindingTimeDist;
	int lysisBlocks;
	int unbindsPerBlock;
	int maxLysesPerBlock;

} Parameters;

typedef struct {
	// input files
	char* dataFileDir;
	char* unbindingTimeFile;
	char* lysisTimeFile;
	char* totalLysesFile;

	// output files
	char* boundMolecules;
	char* unboundMolecules;
	char* fibers;
} Files;

typedef struct {

	double* unbindingTimeDistribution;
	double **lysisTime;
	int* lysesPerBlock;

} InputData;

typedef struct {
	int length;
	int width;
	int total;
} RankGrid;

typedef struct {
	int nodeHeightTotal;
	int nodeLengthTotal;
	int nodeWidthTotal;
	int ghostRows;
} TotalGridSize;

void initializeParameters(Parameters* parameters);
void initializeFiles(Files* files);
void initializeRankGrid(RankGrid* rankGrid);
void initializeTotalGridSize(TotalGridSize* totalGridSize);
void allocateArrays(InputData* inputData, Parameters* parameters);
void setRandomNumberGenerator(int rank, int size);