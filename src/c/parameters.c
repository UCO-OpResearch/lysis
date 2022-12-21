
#include "all.h"

void initializeParameters(Parameters* parameters){

	/******************************************************************************
	** Physical Parameters
	** i.e., those parameters that are measured from the real world
	******************************************************************************/
	// The tPA binding rate. Units of inverse (micromolar*sec)
	parameters->bindingRate = 1.0e-2; // kon

	// Pore size (distance betwen nodes), measured in centimeters
	parameters->poreSize = 1.0135e-4; // delx

	// Diffusion coefficient, measured in cm^2/s
	parameters->diffusionCoeff = 5.0e-7; // Diff

	// Concentration of binding sites in micromolar
	parameters->bindingSites = 4.27e+2; // bs

	/*
	* Distance from the start of one fiber to the next, in microns
	* because distance between nodes is 1.0135 micron
	* and diameter of 1 fiber is 0.0727 micron
	*/
	parameters->gridSymmetryDistance = 1.0862; // dist

	/* The total number of tPA molecules:
	*      43074 is Colin's [tPA]=0.6 nM
	*      86148 is Colin's [tPA]=1.2 nM
	*/
	parameters->totalMolecules = 43074; //43074; //5000; // M (43074)

	// The probability of moving. Make sure it is small enough that we've converged.
	parameters->movingProbability = 0.5; // q (0.2)

	/******************************************************************************
	** Model Time Parameters
	** i.e., those parameters selected to define
	**       the running time and granularity of the model
	******************************************************************************/
	// Total running time for model in seconds
	parameters->totalTime = 1800; //30 * 60; //10 * 60; // tf

	// The length of one timestep, in seconds
	parameters->timeStep = parameters->movingProbability * pow(parameters->poreSize, 2) / (12 * parameters->diffusionCoeff); // tstep

	// The total number of timesteps
	//This may need a ceiling or floor function to return correct int
	parameters->totalTimeSteps = parameters->totalTime / parameters->timeStep;

	/******************************************************************************
	** Data Parameters
	**
	******************************************************************************/

	// Data size parameters
	parameters->unbindingTimeDist = 101;
	parameters->lysisBlocks = 100;
	parameters->unbindsPerBlock = 500;
	parameters->maxLysesPerBlock = 283;

}

void initializeFiles(Files* files) {

	/*
	* lysismat_PLG135_Q2.dat is a matrix with column corresponding to
	* bin number (1-100) and with entries equal to the lysis times obtained
	* in that bin.
	* An entry of 6000 means lysis didn't happen.
	* i.e. the first column, lysismat(:,1), gives the lysis times
	* for the first 100 tPA leaving times.
	*/

	/*
	* lenlysisvect_PLG135_Q2.dat saves the first row entry in each column of
	* lysismat_PLG135_Q2.dat where lysis did not occur,
	* i.e., the first entry there's a 6000
	* i.e., out of the unbinds in the nth percentile (with respect to time),
	*/

	// Data file names
	char dataFileDir[] = "data/";
	char unbindingTimeFile[] = "data/tsectPAPLG135_Q2.dat";
	char lysisTimeFile[] = "data/lysismat_PLG135_Q2.dat";
	char totalLysesFile[] = "data/lenlysisvect_PLG135_Q2.dat";

	// Dynamically allocate
	files->dataFileDir = (char*)malloc(sizeof(dataFileDir) + 1);
	files->unbindingTimeFile = (char*)malloc(sizeof(unbindingTimeFile) + 1);
	files->lysisTimeFile = (char*)malloc(sizeof(lysisTimeFile) + 1);
	files->totalLysesFile = (char*)malloc(sizeof(totalLysesFile) + 1);

	strcpy(files->dataFileDir, dataFileDir);
	strcpy(files->unbindingTimeFile, unbindingTimeFile);
	strcpy(files->lysisTimeFile, lysisTimeFile);
	strcpy(files->totalLysesFile, totalLysesFile);

	// Output file names
	char boundMolecules[] = "output/boundMolecules.dat";
	char unboundMolecules[] = "output/unboundMolecules.dat";
	char fibers[] = "output/fibers.dat";

	// Dynamically allocate
	files->boundMolecules = (char*)malloc(sizeof(boundMolecules) + 1);
	files->unboundMolecules = (char*)malloc(sizeof(unboundMolecules) + 1);
	files->fibers = (char*)malloc(sizeof(fibers) + 1);

	strcpy(files->boundMolecules, boundMolecules);
	strcpy(files->unboundMolecules, unboundMolecules);
	strcpy(files->fibers, fibers);
}

void initializeRankGrid(RankGrid* rankGrid) {
	rankGrid->length = 2;
	rankGrid->width = 7;
	rankGrid->total = rankGrid->length * rankGrid->width; // The total needs to match the number of mpi ranks.
}

// These are total nodes across the entire grid, all ranks included.
// These totals are spread across the mpi ranks.
void initializeTotalGridSize(TotalGridSize* totalGridSize){
	totalGridSize->nodeHeightTotal = 1; // should be 1 unless 3d grid is implemented
	totalGridSize->nodeLengthTotal = 151; // Do not include ghost rows.
	totalGridSize->nodeWidthTotal = 115;
	totalGridSize->ghostRows = 2; // Only applies to bottom boundary ranks. This is added to the nodelength of the bottom boundary ranks.
}

void allocateArrays(InputData* inputData, Parameters* parameters) {

	inputData->unbindingTimeDistribution = (double*)malloc(sizeof(double) * parameters->unbindingTimeDist); // tsec1;

	inputData->lysisTime = (double**)malloc(sizeof(double*) * parameters->maxLysesPerBlock); // lysismat

	for (int iii = 0; iii < parameters->maxLysesPerBlock; iii++) {
		inputData->lysisTime[iii] = (double*)malloc(sizeof(double) * parameters->lysisBlocks);
	}

	inputData->lysesPerBlock = (int*)malloc(sizeof(int) * parameters->lysisBlocks);
}

void setRandomNumberGenerator(int rank, int size) {
	
	// old code
	// Seed for the random number generator
	/*UINT_LEAST32_T seed = 912309035; // seed
	UINT_LEAST32_T state[4] = { 129281, 362436069, 123456789, seed }; // state
	set_kiss32_(state);*/

	if(HEAD_RANK == rank){
		UINT_LEAST32_T** states;
		states = (UINT_LEAST32_T**)malloc(sizeof(UINT_LEAST32_T*) * size);
		for(int iii = 0; iii < size; iii++){
			states[iii] = (UINT_LEAST32_T*)malloc(sizeof(UINT_LEAST32_T) * 4);
		}
		
		//will need to rewrite this when we start using different states for each node
		for(int node = 0; node < size; node++){
			states[node][0] = 129281;
			states[node][1] = 362436069;
			states[node][2] = 123456789;
			states[node][3] = 912309035;
		}

		for(int iii = 1; iii < size; iii++){
			MPI_Send(states[iii], 4, MPI_UINT32_T, iii, 0, MPI_COMM_WORLD);
		}

		set_kiss32_(states[0]);
	} else{
		UINT_LEAST32_T receive[4];
		MPI_Recv(&receive, 4, MPI_UINT32_T, HEAD_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("receive: %i %i %i %i\n", receive[0], receive[1], receive[2], receive[3]);
		//fflush(stdout);
		set_kiss32_(receive);
	}
	//printf("random number: %i %i %i %i\n", urcw1_(), urcw1_(), urcw1_(), urcw1_());
	//fflush(stdout);
	//MPI_Barrier(MPI_COMM_WORLD);
	//MPI_Abort(MPI_COMM_WORLD, 0);
}

