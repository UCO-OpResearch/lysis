#pragma once

#include <math.h>//for UINT_LEAST32_T

 const struct Parameters{

	/******************************************************************************
	** Physical Parameters
	** i.e., those parameters that are measured from the real world
	******************************************************************************/
	// The tPA binding rate. Units of inverse (micromolar*sec)
	double bindingRate; // kon

									   // Pore size (distance betwen nodes), measured in centimeters
	double poreSize; // delx

									   // Diffusion coefficient, measured in cm^2/s
	double diffusionCoeff; // Diff

										  // Concentration of binding sites in micromolar
	double bindingSites; // bs

										 /*
										 * Distance from the start of one fiber to the next, in microns
										 * because distance between nodes is 1.0135 micron
										 * and diameter of 1 fiber is 0.0727 micron
										 */
	double gridSymmetryDistance; // dist

												/* The total number of tPA molecules:
												*      43074 is Colin's [tPA]=0.6 nM
												*      86148 is Colin's [tPA]=1.2 nM
												*/
	int totalMolecules; // M (43074)

									// The probability of moving. Make sure it is small enough that we've converged.
	double movingProbability; // q (0.2)

										  /******************************************************************************
										  ** Model Time Parameters
										  ** i.e., those parameters selected to define
										  **       the running time and granularity of the model
										  ******************************************************************************/
										  // Total running time for model in seconds
	int totalTime; // tf

								   // The length of one timestep, in seconds
	double timeStep; // tstep

																						  // The total number of timesteps
																						  //This may need a ceiling or floor function to return correct int
	int totalTimeSteps;

	/******************************************************************************
	** Data Parameters
	**
	******************************************************************************/

	// Data size parameters
	int unbindingTimeDist;
	int lysisBlocks;
	int unbindsPerBlock;
	int maxLysesPerBlock;

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
	/*std::string dataFileDir = "data/";
	std::string unbindingTimeFile = dataFileDir + "tsectPAPLG135_Q2.dat";
	std::string lysisTimeFile = dataFileDir + "lysismat_PLG135_Q2.dat";
	std::string totalLysesFile = dataFileDir + "lenlysisvect_PLG135_Q2.dat";*/

 } parameters = {
	 .bindingRate = 1.0e-2,
	 .poreSize = 1.0135e-4,
	 .diffusionCoeff = 5.0e-7,
	 .bindingSites = 4.27e+2,
	 .gridSymmetryDistance = 1.0862,
	 .totalMolecules = 500,
	 .movingProbability = 0.5,
	 .totalTime = 10 * 60,
	 .timeStep = movingProbability * pow(poreSize, 2) / (12 * diffusionCoeff),
	 .totalTimeSteps = .totalTime / .timeStep,
	 .unbindingTimeDist = 101,
	 .lysisBlocks = 100,
	 .unbindsPerBlock = 500,
	 .maxLysesPerBlock = 283
 };

typedef struct {

	// The following variables will hold the microscale data
	// which will be read in from the named data files

	double* UnbindingTimeDistribution;
	double **lysisTime;
	int* lysesPerBlock;

} ParameterArrays;

const struct ProcessGrid{
	const int length;
	const int width;
	//const int total;
} processGrid = {
	.length = 4,
	.width = 4
	//.total = .length *.width
};

const struct NodeGridSize {
	const int height;
	const int length;
	const int width;
	int ghostRows;
} nodeGridSize = {
	.height = 1,
	.length = 100,
	.width = 100,
	.ghostRows = 2
};