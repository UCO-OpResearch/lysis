
#include "all.h"

#ifdef BACKTRACE
void backTrace(int signal) {
	void *array[20];
	size_t size;
	size = backtrace(array, 20);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(1);
}
#endif

int main(int argc, char *argv[]) {

	#ifdef BACKTRACE
	signal(SIGSEGV, backTrace);
	#endif

	//************************************************//
	//************************************************//
	// the usual MPI startup routines

	MPI_Init(&argc, &argv);
	int mpiSize, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	RankGrid rankGrid;
	initializeRankGrid(&rankGrid);

	//check that MPI size matches parameter size
	if (mpiSize != rankGrid.total) {
		if (rank == 0) {
			printf("MPI size does not match expected rankGrid.total! MPI Size = %i rankGrid.total = %i\n", mpiSize, rankGrid.total );
		}
		MPI_Abort(MPI_COMM_WORLD, 0);
		return 1;
	}

    char hostname[256];
	gethostname(hostname, sizeof(hostname));
	printf("Hostname: %s Rank: %i\n", hostname, rank);

	//************************************************//
	//************************************************//
	// Setup the simulation

	#ifdef NEWTIMELIST
	if(HEAD_RANK == rank){
		printf("Using new DynArray.\n");
	}
	#endif

	NodeGrid* grid = (NodeGrid*)malloc(sizeof(NodeGrid));
	grid->MPIRank = rank; // needed for setting grid boundaries and many other things
	setGridBoundaries(grid, &rankGrid);

	TotalGridSize totalGridSize;
	initializeTotalGridSize(&totalGridSize);

	Parameters parameters;
	initializeParameters(&parameters);

	Files files;
	initializeFiles(&files);

	//create experiment data arrays
	InputData inputData;
	allocateArrays(&inputData, &parameters);

	//load data from files
	if (0 == rank) {
		if (!loadInputData(&files, &inputData, &parameters)) {
			printf("loadInputData() failure\n");
			fflush(stdout);
			MPI_Abort(MPI_COMM_WORLD, 0);
			return 1;
		}
	}

	transferInputData(&inputData, &parameters);

	setRandomNumberGenerator(rank, mpiSize);

	//create grid
	initializeNodeGrid(grid, &totalGridSize, &parameters, &rankGrid);

	// printf("My rank = %i length = %i lengthWithGhostRows = %i width = %i, height = %i\n", rank, grid->nodeLength, grid->nodeLengthWithGhostRows, grid->nodeWidth, grid->nodeHeight);
	// printf("My rank = %i length %i width = %i\n", rank, grid->nodeLength,  grid->nodeWidth);
	// MPI_Finalize();
	// return 0;

	// head rank will randomly distribute the molecules. Then it will hand them out to the other bottom ranks.
	distributeMolecules(grid, &rankGrid, parameters.totalMolecules);

	//************************************************//
	//************************************************//
	// run simulation

	double startTime;
	if (0 == rank) {
		startTime = MPI_Wtime();
	}

	simulationLoop(grid, &parameters, &inputData, &rankGrid, &files);

	if (0 == rank) {
		double endTime;
		endTime = MPI_Wtime();
		printf("Run time: %lf minutes\n", (endTime - startTime)/60);
	}

	//************************************************//
	//************************************************//
	
	free(grid);
	MPI_Finalize();
	return 0;
}